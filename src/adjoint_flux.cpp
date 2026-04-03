#include "openmc/adjoint_flux.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"

namespace openmc {

void AdjointFlux::compute_from_files(const std::string& transfer_function_file,
  const std::string& fission_matrix_file, const std::string& output_file)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "共轭通量计算" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  // 1. 读取传递函数数据
  std::cout << "\n读取传递函数: " << transfer_function_file << std::endl;

  // 检查传递函数文件是否存在
  if (!file_exists(transfer_function_file)) {
    fatal_error("传递函数文件不存在: " + transfer_function_file +
                "\n请确保已计算并保存传递函数。");
  }

  hid_t tf_file = file_open(transfer_function_file, 'r');

  // 读取网格参数
  read_dataset(tf_file, "origin", origin_);
  read_dataset(tf_file, "shape", shape_);
  read_attribute(tf_file, "pitch", pitch_);

  int tf_n_groups = 1;
  int tf_n_families = 1;
  vector<double> tf_energy_edges;
  if (object_exists(tf_file, "n_groups")) {
    read_dataset(tf_file, "n_groups", tf_n_groups);
  }
  if (object_exists(tf_file, "n_families")) {
    read_dataset(tf_file, "n_families", tf_n_families);
  }
  if (object_exists(tf_file, "energy_edges")) {
    read_dataset(tf_file, "energy_edges", tf_energy_edges);
  }
  if (!tf_energy_edges.empty() &&
      static_cast<int>(tf_energy_edges.size()) - 1 != tf_n_groups) {
    fatal_error("AdjointFlux: n_groups does not match energy_edges length in "
                "transfer function file.");
  }

  n_cells_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  std::cout << "  网格: " << shape_[0] << "x" << shape_[1] << "x" << shape_[2]
            << "=" << n_cells_ << "单元, 间距: " << pitch_ << " cm"
            << std::endl;

  // 读取源单元索引列表
  vector<int> source_cell_indices;
  read_dataset(tf_file, "source_cell_indices", source_cell_indices);

  // 读取传递函数（稀疏格式）
  std::unordered_map<int, std::unordered_map<int, vector<double>>>
    transfer_functions;

  hid_t tf_group = H5Gopen(tf_file, "transfer_functions", H5P_DEFAULT);

  size_t total_tf_entries = 0;
  for (int i_source : source_cell_indices) {
    std::string cell_name = "source_cell_" + std::to_string(i_source);

    if (object_exists(tf_group, cell_name.c_str())) {
      hid_t cell_group = H5Gopen(tf_group, cell_name.c_str(), H5P_DEFAULT);

      vector<int> indices;
      vector<double> values;
      read_dataset(cell_group, "indices", indices);
      read_dataset(cell_group, "values", values);

      size_t entry_count = indices.size();
      size_t entries_per_response =
        static_cast<size_t>(tf_n_families) * tf_n_groups;
      size_t expected = entries_per_response * entry_count;
      if (tf_n_families == 1 && tf_n_groups == 1 &&
          values.size() == entry_count) {
        expected = entry_count;
      }

      if (values.size() != expected) {
        fatal_error("AdjointFlux: transfer function values dataset has an "
                    "unexpected length.");
      }

      // 构建稀疏 map（每个响应单元存 n_families * n_groups 个值）
      for (size_t k = 0; k < indices.size(); ++k) {
        vector<double> fg_values(
          static_cast<size_t>(entries_per_response), 0.0);
        if (entries_per_response == 1 && values.size() == entry_count) {
          fg_values[0] = values[k];
        } else {
          for (size_t v = 0; v < entries_per_response; ++v) {
            size_t offset = k * entries_per_response + v;
            fg_values[v] = values[offset];
          }
        }
        transfer_functions[i_source][indices[k]] = std::move(fg_values);
      }

      total_tf_entries += indices.size();
      H5Gclose(cell_group);
    }
  }

  H5Gclose(tf_group);
  file_close(tf_file);

  std::cout << "  源单元: " << source_cell_indices.size()
            << ", 传递函数条目: " << total_tf_entries << std::endl;

  // 2. 读取伴随源数据
  std::cout << "\n读取伴随源: " << fission_matrix_file << std::endl;

  // 检查裂变矩阵文件是否存在
  if (!file_exists(fission_matrix_file)) {
    fatal_error("裂变矩阵文件不存在: " + fission_matrix_file +
                "\n请确保已计算并保存裂变矩阵。");
  }

  hid_t fm_file = file_open(fission_matrix_file, 'r');

  // 检查是否存在伴随源
  if (!object_exists(fm_file, "adjoint_source")) {
    file_close(fm_file);
    fatal_error(
      "裂变矩阵文件中未找到伴随源。请先运行 compute_adjoint_source。");
  }

  vector<double> adjoint_source;
  read_dataset(fm_file, "adjoint_source", adjoint_source);

  double keff_reference = 0.0;
  if (attribute_exists(fm_file, "keff_reference")) {
    read_attribute(fm_file, "keff_reference", keff_reference);
  }

  int adjoint_iterations = 0;
  if (attribute_exists(fm_file, "adjoint_iterations")) {
    read_attribute(fm_file, "adjoint_iterations", adjoint_iterations);
  }

  file_close(fm_file);

  std::cout << "  伴随源单元: " << adjoint_source.size();
  if (keff_reference > 0.0) {
    std::cout << ", keff=" << std::fixed << std::setprecision(5)
              << keff_reference;
  }
  std::cout << std::endl;

  // 验证网格一致性
  if (adjoint_source.size() != n_cells_) {
    fatal_error("网格不匹配: 传递函数有 " + std::to_string(n_cells_) +
                " 个单元，但伴随源有 " + std::to_string(adjoint_source.size()) +
                " 个单元");
  }

  // 3. 计算共轭通量
  std::cout << "\n计算共轭通量: Φ†(r) = Σ_i T(i→r) × S†(i)" << std::endl;

  compute_from_memory(transfer_functions, adjoint_source, shape_, origin_,
    pitch_, tf_n_groups, tf_n_families, std::move(tf_energy_edges));

  // 4. 写入输出文件
  write_to_file(output_file);

  std::cout << std::string(70, '=') << std::endl;
}

void AdjointFlux::compute_from_memory(
  const std::unordered_map<int, std::unordered_map<int, vector<double>>>&
    transfer_functions,
  const vector<double>& adjoint_source, const std::array<int, 3>& shape,
  const std::array<double, 3>& origin, double pitch, int n_groups,
  int n_families, vector<double> energy_edges)
{
  shape_ = shape;
  origin_ = origin;
  pitch_ = pitch;
  n_cells_ = static_cast<size_t>(shape[0]) * shape[1] * shape[2];
  n_groups_ = n_groups > 0 ? n_groups : 1;
  n_families_ = n_families > 1 ? n_families : 1;
  has_family_data_ = (n_families_ > 1);
  energy_edges_ = std::move(energy_edges);
  if (n_groups_ > 1) {
    if (energy_edges_.size() != static_cast<size_t>(n_groups_ + 1)) {
      fatal_error("AdjointFlux: energy_edges must have n_groups + 1 entries.");
    }
    if (!std::is_sorted(energy_edges_.begin(), energy_edges_.end())) {
      fatal_error("AdjointFlux: energy_edges must be sorted ascending.");
    }
  } else {
    n_groups_ = 1;
    energy_edges_.clear();
  }

  group_total_flux_.assign(static_cast<size_t>(n_groups_), 0.0);

  // 清空之前的结果
  adjoint_flux_sparse_.clear();
  family_adjoint_flux_.clear();
  if (has_family_data_) {
    family_adjoint_flux_.resize(n_families_);
  }

  // 计算 Φ†(r) = Σ_i T(i -> r) × S†(i)
  // 当 has_family_data_ 时，同时计算 family-resolved 和 total 卷积
  for (const auto& [i_source, response_map] : transfer_functions) {
    if (i_source < 0 || i_source >= static_cast<int>(adjoint_source.size())) {
      continue;
    }

    double importance = adjoint_source[i_source];
    if (importance <= 0.0) {
      continue;
    }

    for (const auto& [j_response, T_values] : response_map) {
      auto& total_flux = adjoint_flux_sparse_[j_response];
      if (total_flux.empty()) {
        total_flux = make_zero_group_vector();
      }

      if (has_family_data_) {
        // Family-resolved: T_values has n_families * n_groups entries
        for (int f = 0; f < n_families_; ++f) {
          auto& family_flux = family_adjoint_flux_[f][j_response];
          if (family_flux.empty()) {
            family_flux = make_zero_group_vector();
          }
          for (int g = 0; g < n_groups_; ++g) {
            size_t idx = static_cast<size_t>(f) * n_groups_ + g;
            double contrib = T_values[idx] * importance;
            family_flux[g] += contrib;
            total_flux[g] += contrib;
          }
        }
      } else {
        // No family data: T_values has n_groups entries (backward compatible)
        for (int g = 0; g < n_groups_; ++g) {
          total_flux[g] += T_values[g] * importance;
        }
      }
    }
  }

  // 计算统计信息
  max_flux_ = 0.0;
  total_flux_ = 0.0;
  nonzero_cells_ = 0;

  for (const auto& [j, flux_vector] : adjoint_flux_sparse_) {
    double cell_total = 0.0;
    for (int g = 0; g < n_groups_; ++g) {
      group_total_flux_[g] += flux_vector[g];
      cell_total += flux_vector[g];
    }
    if (cell_total > 0.0) {
      ++nonzero_cells_;
    }
    max_flux_ = std::max(max_flux_, cell_total);
    total_flux_ += cell_total;
  }

  // 输出统计信息
  double density = 100.0 * nonzero_cells_ / n_cells_;
  std::cout << "\n共轭通量统计:" << std::endl;
  std::cout << "  非零单元: " << nonzero_cells_ << "/" << n_cells_ << " ("
            << std::fixed << std::setprecision(1) << density << "%)"
            << ", 最大值: " << std::scientific << std::setprecision(3)
            << max_flux_ << ", 总通量: " << total_flux_ << std::endl;
}

vector<double> AdjointFlux::get_adjoint_flux_dense() const
{
  vector<double> dense(n_cells_, 0.0);

  for (const auto& [j, flux_vector] : adjoint_flux_sparse_) {
    if (j >= 0 && j < static_cast<int>(n_cells_)) {
      dense[j] = std::accumulate(flux_vector.begin(), flux_vector.end(), 0.0);
    }
  }

  return dense;
}

vector<double> AdjointFlux::get_adjoint_flux_group_dense() const
{
  vector<double> dense(static_cast<size_t>(n_groups_) * n_cells_, 0.0);

  for (const auto& [j, flux_vector] : adjoint_flux_sparse_) {
    if (j < 0 || j >= static_cast<int>(n_cells_))
      continue;

    for (int g = 0; g < n_groups_; ++g) {
      size_t idx = static_cast<size_t>(g) * n_cells_ + j;
      dense[idx] = flux_vector[g];
    }
  }

  return dense;
}

double AdjointFlux::get_max_value() const
{
  return max_flux_;
}

double AdjointFlux::get_total_flux() const
{
  return total_flux_;
}

void AdjointFlux::write_to_file(const std::string& filename)
{
  std::cout << "\n写入共轭通量: " << filename << std::endl;

  hid_t file_id = file_open(filename, 'w');

  // 写入文件属性
  write_attribute(file_id, "filetype", "adjoint_flux");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "description",
    "Adjoint flux computed from transfer function and adjoint source");
  write_attribute(file_id, "storage_format", "sparse");

  // 写入网格信息
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "grid_shape", shape_);
  std::array<double, 3> grid_pitch_array {pitch_, pitch_, pitch_};
  write_dataset(file_id, "grid_pitch", grid_pitch_array);
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells_));
  write_dataset(file_id, "n_groups", n_groups_);
  if (!energy_edges_.empty()) {
    write_dataset(file_id, "energy_edges", energy_edges_);
  }

  // 写入统计信息
  write_attribute(file_id, "nonzero_cells", static_cast<int>(nonzero_cells_));
  write_attribute(file_id, "max_flux", max_flux_);
  write_attribute(file_id, "total_flux", total_flux_);
  write_attribute(
    file_id, "density_percent", 100.0 * nonzero_cells_ / n_cells_);

  // 写入稀疏格式数据
  vector<int> indices;
  vector<double> total_values;
  vector<double> group_values;

  indices.reserve(adjoint_flux_sparse_.size());
  total_values.reserve(adjoint_flux_sparse_.size());
  group_values.reserve(
    adjoint_flux_sparse_.size() * static_cast<size_t>(n_groups_));

  for (const auto& [j, flux_vector] : adjoint_flux_sparse_) {
    indices.push_back(j);
    double total = std::accumulate(flux_vector.begin(), flux_vector.end(), 0.0);
    total_values.push_back(total);
    group_values.insert(
      group_values.end(), flux_vector.begin(), flux_vector.end());
  }

  // 按索引排序
  vector<size_t> sort_indices(indices.size());
  std::iota(sort_indices.begin(), sort_indices.end(), 0);
  std::sort(sort_indices.begin(), sort_indices.end(),
    [&indices](size_t a, size_t b) { return indices[a] < indices[b]; });

  vector<int> sorted_indices(indices.size());
  vector<double> sorted_values(total_values.size());
  for (size_t i = 0; i < sort_indices.size(); ++i) {
    sorted_indices[i] = indices[sort_indices[i]];
    sorted_values[i] = total_values[sort_indices[i]];
  }

  write_dataset(file_id, "cell_indices", sorted_indices);
  write_dataset(file_id, "flux_mean", sorted_values);

  // 写入分群数据
  vector<double> sorted_group_values(group_values.size());
  for (size_t i = 0; i < sort_indices.size(); ++i) {
    size_t src_offset = sort_indices[i] * static_cast<size_t>(n_groups_);
    size_t dst_offset = i * static_cast<size_t>(n_groups_);
    std::copy_n(group_values.begin() + src_offset, n_groups_,
      sorted_group_values.begin() + dst_offset);
  }
  write_dataset(file_id, "flux_group_mean", sorted_group_values);

  // 写入稠密格式
  auto dense_flux = get_adjoint_flux_dense();
  write_dataset(file_id, "adjoint_flux_dense", dense_flux);
  auto dense_group_flux = get_adjoint_flux_group_dense();
  write_dataset(file_id, "adjoint_flux_group_dense", dense_group_flux);

  if (!group_total_flux_.empty()) {
    write_dataset(file_id, "group_total_flux", group_total_flux_);
  }

  // ========== family_resolved/ 组 ==========
  if (has_family_data_ && !family_adjoint_flux_.empty()) {
    hid_t fr_group = H5Gcreate(
      file_id, "family_resolved", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write_attribute(fr_group, "n_families", n_families_);
    write_attribute(fr_group, "description",
      "Family-resolved adjoint flux: phi_dag_f(r) = sum_i T_f(i->r) * I*(i)");
    write_attribute(
      fr_group, "family_order", "prompt, delayed_1, ..., delayed_8");

    const char* family_names[] = {"prompt", "delayed_1", "delayed_2",
      "delayed_3", "delayed_4", "delayed_5", "delayed_6", "delayed_7",
      "delayed_8"};

    for (int f = 0; f < n_families_ && f < 9; ++f) {
      const auto& family_data = family_adjoint_flux_[f];
      if (family_data.empty())
        continue;

      hid_t fam_group = H5Gcreate(
        fr_group, family_names[f], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Collect data
      vector<int> f_indices;
      vector<double> f_mean;
      vector<double> f_group_vals;

      for (const auto& [j, fv] : family_data) {
        f_indices.push_back(j);
        double total = std::accumulate(fv.begin(), fv.end(), 0.0);
        f_mean.push_back(total);
        f_group_vals.insert(f_group_vals.end(), fv.begin(), fv.end());
      }

      // Sort by cell index
      vector<size_t> si(f_indices.size());
      std::iota(si.begin(), si.end(), 0);
      std::sort(si.begin(), si.end(), [&f_indices](size_t a, size_t b) {
        return f_indices[a] < f_indices[b];
      });

      vector<int> s_idx(f_indices.size());
      vector<double> s_mean(f_mean.size());
      vector<double> s_grp(f_group_vals.size());
      for (size_t i = 0; i < si.size(); ++i) {
        s_idx[i] = f_indices[si[i]];
        s_mean[i] = f_mean[si[i]];
        size_t src = si[i] * static_cast<size_t>(n_groups_);
        size_t dst = i * static_cast<size_t>(n_groups_);
        std::copy_n(f_group_vals.begin() + src, n_groups_, s_grp.begin() + dst);
      }

      write_dataset(fam_group, "cell_indices", s_idx);
      write_dataset(fam_group, "flux_mean", s_mean);
      if (n_groups_ > 1) {
        write_dataset(fam_group, "flux_group_mean", s_grp);
      }
      write_attribute(fam_group, "n_cells", static_cast<int>(f_indices.size()));

      H5Gclose(fam_group);
    }

    H5Gclose(fr_group);

    std::cout << "  Family-resolved: " << n_families_ << " families written"
              << std::endl;
  }

  // ========== 语义元数据 (Phase 2) ==========
  // 记录此文件中 "共轭通量" 的物理含义
  hid_t sem_group = H5Gcreate(
    file_id, "semantic_metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  write_attribute(sem_group, "physical_quantity",
    "Response-weighted importance field, NOT true adjoint flux");
  write_attribute(sem_group, "definition",
    "phi_dag(r) = sum_i T(i->r) * I_star(i), where T is the transfer "
    "function and I_star is the fission matrix eigenvector (scalar per cell)");
  write_attribute(sem_group, "group_meaning",
    "Group data phi_dag_g(r) reflects collision-energy-resolved importance "
    "at response position r. The group index g refers to the energy of the "
    "neutron causing fission at r, not the birth energy of fission neutrons.");
  write_attribute(sem_group, "adjoint_source",
    "I_star(i) from fission matrix eigenvalue problem. Scalar (no energy "
    "dependence). This is the main limitation for beta_eff accuracy.");
  write_attribute(sem_group, "limitations",
    "1) I_star is scalar, not energy-resolved => cannot distinguish prompt/"
    "delayed importance at source level. 2) phi_dag_g is response-energy "
    "(collision), not birth-energy => chi weighting is approximate.");
  write_attribute(sem_group, "usage_in_beta_eff",
    "Method A: weight by chi (birth spectrum). "
    "Method B: use scalar I_star directly. "
    "Method D: weight by nu fractions (production decomposition).");
  H5Gclose(sem_group);

  file_close(file_id);

  std::cout << "  稀疏条目: " << adjoint_flux_sparse_.size() << ", 估计大小: "
            << (adjoint_flux_sparse_.size() * 12 + n_cells_ * 8) / (1024 * 1024)
            << " MB" << std::endl;
}

std::vector<double> AdjointFlux::make_zero_group_vector() const
{
  return std::vector<double>(static_cast<size_t>(n_groups_), 0.0);
}

} // namespace openmc
