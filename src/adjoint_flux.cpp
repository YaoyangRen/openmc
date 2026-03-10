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
  vector<double> tf_energy_edges;
  if (object_exists(tf_file, "n_groups")) {
    read_dataset(tf_file, "n_groups", tf_n_groups);
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
      size_t expected = static_cast<size_t>(tf_n_groups) * entry_count;
      if (tf_n_groups == 1 && values.size() == entry_count) {
        expected = entry_count;
      }

      if (values.size() != expected) {
        fatal_error("AdjointFlux: transfer function values dataset has an "
                    "unexpected length.");
      }

      // 构建稀疏 map（按能群存储）
      for (size_t k = 0; k < indices.size(); ++k) {
        vector<double> group_values(static_cast<size_t>(tf_n_groups), 0.0);
        if (tf_n_groups == 1 && values.size() == entry_count) {
          group_values[0] = values[k];
        } else {
          for (int g = 0; g < tf_n_groups; ++g) {
            size_t offset = static_cast<size_t>(k) * tf_n_groups + g;
            group_values[g] = values[offset];
          }
        }
        transfer_functions[i_source][indices[k]] = std::move(group_values);
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
    pitch_, tf_n_groups, std::move(tf_energy_edges));

  // 4. 写入输出文件
  write_to_file(output_file);

  std::cout << std::string(70, '=') << std::endl;
}

void AdjointFlux::compute_from_memory(
  const std::unordered_map<int, std::unordered_map<int, vector<double>>>&
    transfer_functions,
  const vector<double>& adjoint_source, const std::array<int, 3>& shape,
  const std::array<double, 3>& origin, double pitch, int n_groups,
  vector<double> energy_edges)
{
  shape_ = shape;
  origin_ = origin;
  pitch_ = pitch;
  n_cells_ = static_cast<size_t>(shape[0]) * shape[1] * shape[2];
  n_groups_ = n_groups > 0 ? n_groups : 1;
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

  // 计算 Φ†(r) = Σ_i T(i -> r) × S†(i)
  for (const auto& [i_source, response_map] : transfer_functions) {
    if (i_source < 0 || i_source >= static_cast<int>(adjoint_source.size())) {
      continue;
    }

    double importance = adjoint_source[i_source];
    if (importance <= 0.0) {
      continue;
    }

    for (const auto& [j_response, T_values] : response_map) {
      auto& flux_vector = adjoint_flux_sparse_[j_response];
      if (flux_vector.empty()) {
        flux_vector = make_zero_group_vector();
      }
      for (int g = 0; g < n_groups_; ++g) {
        flux_vector[g] += T_values[g] * importance;
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
