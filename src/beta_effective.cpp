#include "openmc/beta_effective.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/material_nuclear_data.h"
#include "openmc/nuclide.h"
#include "openmc/position.h"
#include "openmc/reaction_product.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <vector>

namespace openmc {

//------------------------------------------------------------------------------
// 核心驱动函数：根据网格通量/共轭通量文件完成 β_eff 计算。
// 算法遵循扰动理论定义：
//   β_i,eff = ∫ φ*(r,E) χ_d,i(E) ν_d,i(E,r) Σ_f(r,E) φ(r,E) dV dE / 分母
// 分母使用瞬发量。实现直接采用材料相关核数据路径。
//------------------------------------------------------------------------------
void BetaEffective::compute_from_files(const std::string& flux_file,
  const std::string& adjoint_flux_file, const std::string& output_file)
{
  std::cout << std::string(70, '=') << std::endl;

  // 1. 读取正向通量
  std::cout << "[1/4] 读取正向通量: " << flux_file << std::endl;
  if (!file_exists(flux_file)) {
    fatal_error("通量文件不存在: " + flux_file + "\n请先运行正向计算。");
  }

  std::unordered_map<int, double> flux;
  std::array<int, 3> flux_shape {};
  double flux_pitch = 0.0;
  read_flux_data(flux_file, flux, flux_shape, flux_pitch);

  std::cout << "  网格: " << flux_shape[0] << "x" << flux_shape[1] << "x"
            << flux_shape[2] << ", 间距: " << flux_pitch << " cm"
            << ", 非零单元: " << flux.size() << std::endl;

  // 2. 读取共轭通量 (φ*)
  std::cout << "\n[2/4] 读取共轭通量: " << adjoint_flux_file << std::endl;
  if (!file_exists(adjoint_flux_file)) {
    fatal_error(
      "共轭通量文件不存在: " + adjoint_flux_file + "\n请确保已计算共轭通量。");
  }

  std::unordered_map<int, double> adjoint_flux;
  std::array<int, 3> adjoint_shape {};
  double adjoint_pitch = 0.0;
  read_adjoint_flux_data(
    adjoint_flux_file, adjoint_flux, adjoint_shape, adjoint_pitch);

  std::cout << "  网格: " << adjoint_shape[0] << "x" << adjoint_shape[1] << "x"
            << adjoint_shape[2] << ", 间距: " << adjoint_pitch << " cm"
            << ", 非零单元: " << adjoint_flux.size() << std::endl;

  // 3. 验证网格一致性：β_eff 依赖逐单元的 φ-φ* 内积，因此两张网格必须完全
  //    对齐（尺寸、pitch 以及多群能量边界）。
  if (flux_shape != adjoint_shape ||
      std::abs(flux_pitch - adjoint_pitch) > 1e-6) {
    fatal_error(
      "Grid mismatch between flux and adjoint flux!\n" +
      std::string("  Flux grid: ") + std::to_string(flux_shape[0]) + "x" +
      std::to_string(flux_shape[1]) + "x" + std::to_string(flux_shape[2]) +
      ", pitch=" + std::to_string(flux_pitch) + "\n" +
      std::string("  Adjoint grid: ") + std::to_string(adjoint_shape[0]) + "x" +
      std::to_string(adjoint_shape[1]) + "x" +
      std::to_string(adjoint_shape[2]) +
      ", pitch=" + std::to_string(adjoint_pitch));
  }

  grid_shape_ = flux_shape;
  grid_pitch_ = flux_pitch;

  // 校验多群信息，确保输入网格一致
  validate_group_metadata();
  initialize_flux_spectrum_weights();

  //============================这段代码用于调试=================================
  // 调试输出：基于核数据的材料 β_i（不含通量和共轭通量权重）
  // 用户可以根据材料 id 调整此处参数，用于与 MCNP 的材料 β_i 对比。
  // 这里只打印一次，不影响后续 β_eff 计算。
  int debug_material_id = 1; // TODO: 根据模型实际材料 id 修改
  try {
    auto debug_data = extract_material_nuclear_data(debug_material_id);
    if (debug_data.is_fissionable && !debug_data.sigma_f_groups.empty() &&
        !debug_data.nu_total_groups.empty() &&
        !debug_data.nu_delayed_groups.empty()) {
      const int G = static_cast<int>(debug_data.sigma_f_groups.size());
      const int I = N_DELAYED_GROUPS;

      std::vector<double> num(I, 0.0);
      double den = 0.0;

      for (int g = 0; g < G; ++g) {
        double SigmaF = debug_data.sigma_f_groups[g];
        double nu_tot = debug_data.nu_total_groups[g];
        if (SigmaF <= 0.0)
          continue;

        den += nu_tot * SigmaF;

        for (int i = 0; i < I; ++i) {
          int off = delayed_offset(g, i);
          if (off >= 0 &&
              off < static_cast<int>(debug_data.nu_delayed_groups.size())) {
            double nu_d = debug_data.nu_delayed_groups[off];
            num[i] += nu_d * SigmaF;
          }
        }
      }

      std::cout << "\n[DEBUG] Library-based beta_i for material "
                << debug_data.material_name << " (id=" << debug_material_id
                << ")" << std::endl;
      double beta_sum = 0.0;
      for (int i = 0; i < I; ++i) {
        double beta_i = (den > 0.0) ? num[i] / den : 0.0;
        beta_sum += beta_i;
        std::cout << "  group " << (i + 1)
                  << ": beta_i = " << std::setprecision(8) << beta_i
                  << std::endl;
      }
      std::cout << "  total beta = " << std::setprecision(8) << beta_sum << "\n"
                << std::endl;
    } else {
      std::cout << "\n[DEBUG] Material id " << debug_material_id
                << " is non-fissionable or missing group data; "
                << "skip library beta_i debug.\n";
    }
  } catch (const std::exception& e) {
    std::cout
      << "\n[DEBUG] Failed to compute library-based beta_i for material "
      << debug_material_id << ": " << e.what() << "\n";
  }
  //=============================================================================

  std::fill(beta_i_.begin(), beta_i_.end(), 0.0);
  std::fill(numerators_.begin(), numerators_.end(), 0.0);
  denominator_ = 0.0;

  // 计算网格单元体积
  double volume = flux_pitch * flux_pitch * flux_pitch;

  std::cout << "\n[3/4] 计算 β_eff (材料相关核数据)" << std::endl;
  std::cout << "  采样模式: " << (n_sample_points_ == 1 ? "单点" : "多点")
            << ", 单元体积: " << volume << " cm³" << std::endl;

  // 构建单元-材料映射并提取核数据：通过几何查询找出每个网格实际包含的
  // 材料组合，从而得到空间依赖的 Σ_f、ν_t、ν_d、χ。
  build_cell_material_map(flux);

  // 计算分母 (含 chi_prompt)
  denominator_ = compute_denominator_material(flux, adjoint_flux, volume);

  // 计算不含 chi 的分母作为对比
  double denominator_no_chi =
    compute_denominator_without_chi(flux, adjoint_flux, volume);

  if (denominator_ <= 0.0) {
    fatal_error(
      "Denominator (with chi) is zero or negative! Cannot compute β_eff.");
  }

  std::cout << "\n  分母对比: D(含χ)=" << std::scientific
            << std::setprecision(4) << denominator_
            << ", D(无χ)=" << denominator_no_chi << ", 比值=" << std::fixed
            << std::setprecision(4)
            << (denominator_no_chi > 0 ? denominator_ / denominator_no_chi
                                       : 0.0)
            << std::endl;

  // 计算各缓发群的β_i
  for (int i = 0; i < 8; ++i) {
    numerators_[i] =
      compute_delayed_numerator_material(i, flux, adjoint_flux, volume);
    beta_i_[i] = numerators_[i] / denominator_;
  }
  beta_total_ = std::accumulate(beta_i_.begin(), beta_i_.end(), 0.0);

  // MCNP参考值 (beta_eff)
  const double mcnp_beta[6] = {
    0.00016, 0.00104, 0.00097, 0.00253, 0.00107, 0.00042};
  const double mcnp_beta_total = 0.00619;

  // 输出与MCNP对比结果
  std::cout << "\n  缓发中子先驱核群 β_eff 计算结果与MCNP对比:" << std::endl;
  std::cout << "  " << std::string(60, '-') << std::endl;
  std::cout << "    先驱核群    OpenMC β_eff    MCNP β_eff    相对偏差(%)"
            << std::endl;
  std::cout << "  " << std::string(60, '-') << std::endl;

  for (int i = 0; i < 6; ++i) {
    double rel_err = (mcnp_beta[i] > 0)
                       ? (beta_i_[i] - mcnp_beta[i]) / mcnp_beta[i] * 100.0
                       : 0.0;
    std::cout << "       " << (i + 1) << "        " << std::scientific
              << std::setprecision(5) << beta_i_[i] << "     "
              << std::setprecision(5) << mcnp_beta[i] << "     " << std::fixed
              << std::setprecision(2) << std::setw(7) << rel_err << std::endl;
  }
  // 群7和群8（如有）
  for (int i = 6; i < 8; ++i) {
    if (beta_i_[i] > 1e-10) {
      std::cout << "       " << (i + 1) << "        " << std::scientific
                << std::setprecision(5) << beta_i_[i]
                << "        -             -" << std::endl;
    }
  }

  std::cout << "  " << std::string(60, '-') << std::endl;
  double total_rel_err =
    (mcnp_beta_total > 0)
      ? (beta_total_ - mcnp_beta_total) / mcnp_beta_total * 100.0
      : 0.0;
  std::cout << "      总计      " << std::scientific << std::setprecision(5)
            << beta_total_ << "     " << std::setprecision(5) << mcnp_beta_total
            << "     " << std::fixed << std::setprecision(2) << std::setw(7)
            << total_rel_err << std::endl;
  std::cout << "  " << std::string(60, '-') << std::endl;

  // 7. 输出结果到文件
  std::cout << "\n[4/4] 写入结果: " << output_file << std::endl;
  write_to_file(output_file);

  std::cout << std::string(70, '=') << std::endl;
}

//------------------------------------------------------------------------------

void BetaEffective::read_flux_data(const std::string& filename,
  std::unordered_map<int, double>& flux_map, std::array<int, 3>& shape,
  double& pitch)
{
  hid_t file_id = file_open(filename, 'r');

  flux_group_map_.clear();
  flux_energy_edges_.clear();
  flux_has_group_data_ = false;
  flux_n_groups_ = 1;

  // 读取网格参数
  std::array<int, 3> grid_shape;
  read_dataset(file_id, "grid_shape", grid_shape);
  shape = grid_shape;

  std::array<double, 3> grid_pitch_array;
  read_dataset(file_id, "grid_pitch", grid_pitch_array);
  pitch = grid_pitch_array[0]; // 假设各向同性

  // 读取网格左下角坐标
  if (object_exists(file_id, "grid_lower_left")) {
    read_dataset(file_id, "grid_lower_left", grid_lower_left_);
  }

  if (object_exists(file_id, "n_groups")) {
    read_dataset(file_id, "n_groups", flux_n_groups_);
  }
  if (object_exists(file_id, "energy_edges")) {
    read_dataset(file_id, "energy_edges", flux_energy_edges_);
  }

  // 读取稀疏通量数据
  std::vector<int> cell_indices;
  std::vector<double> flux_mean;
  std::vector<double> flux_group_mean;

  read_dataset(file_id, "cell_indices", cell_indices);
  read_dataset(file_id, "flux_mean", flux_mean);

  if (flux_n_groups_ > 1 && object_exists(file_id, "flux_group_mean")) {
    read_dataset(file_id, "flux_group_mean", flux_group_mean);
    size_t expected = static_cast<size_t>(flux_n_groups_) * cell_indices.size();
    if (flux_group_mean.size() != expected) {
      fatal_error(
        "flux_group_mean size mismatch: expected n_groups * non-zero cells");
    }
    flux_has_group_data_ = true;
  }

  // 构建稀疏 map
  flux_map.clear();
  for (size_t i = 0; i < cell_indices.size(); ++i) {
    flux_map[cell_indices[i]] = flux_mean[i];

    if (flux_has_group_data_) {
      const int cell = cell_indices[i];
      std::vector<double> groups(flux_n_groups_, 0.0);
      size_t offset = static_cast<size_t>(i) * flux_n_groups_;
      for (int g = 0; g < flux_n_groups_; ++g) {
        groups[g] = flux_group_mean[offset + g];
      }
      flux_group_map_[cell] = std::move(groups);
    }
  }

  file_close(file_id);
}

//------------------------------------------------------------------------------

void BetaEffective::read_adjoint_flux_data(const std::string& filename,
  std::unordered_map<int, double>& adjoint_flux_map, std::array<int, 3>& shape,
  double& pitch)
{
  hid_t file_id = file_open(filename, 'r');

  adjoint_group_map_.clear();
  adjoint_energy_edges_.clear();
  adjoint_has_group_data_ = false;
  adjoint_n_groups_ = 1;

  // 读取网格参数
  std::array<int, 3> grid_shape;
  if (object_exists(file_id, "grid_shape")) {
    read_dataset(file_id, "grid_shape", grid_shape);
  } else {
    read_dataset(file_id, "shape", grid_shape);
  }
  shape = grid_shape;

  if (object_exists(file_id, "grid_pitch")) {
    std::array<double, 3> grid_pitch_array;
    read_dataset(file_id, "grid_pitch", grid_pitch_array);
    pitch = grid_pitch_array[0];
  } else {
    read_attribute(file_id, "pitch", pitch);
  }

  if (object_exists(file_id, "n_groups")) {
    read_dataset(file_id, "n_groups", adjoint_n_groups_);
  }
  if (object_exists(file_id, "energy_edges")) {
    read_dataset(file_id, "energy_edges", adjoint_energy_edges_);
  }

  // 读取稀疏共轭通量数据
  std::vector<int> cell_indices;
  std::vector<double> adjoint_flux_values;
  std::vector<double> adjoint_group_values;

  read_dataset(file_id, "cell_indices", cell_indices);
  if (object_exists(file_id, "adjoint_flux_values")) {
    read_dataset(file_id, "adjoint_flux_values", adjoint_flux_values);
  } else {
    read_dataset(file_id, "flux_mean", adjoint_flux_values);
  }

  if (adjoint_n_groups_ > 1) {
    if (object_exists(file_id, "adjoint_flux_group_values")) {
      read_dataset(file_id, "adjoint_flux_group_values", adjoint_group_values);
    } else if (object_exists(file_id, "flux_group_mean")) {
      read_dataset(file_id, "flux_group_mean", adjoint_group_values);
    }
  }

  if (!adjoint_group_values.empty()) {
    size_t expected =
      static_cast<size_t>(adjoint_n_groups_) * cell_indices.size();
    if (adjoint_group_values.size() != expected) {
      fatal_error("adjoint_flux_group_values size mismatch with n_groups");
    }
    adjoint_has_group_data_ = true;
  }

  // 构建稀疏 map
  adjoint_flux_map.clear();
  for (size_t i = 0; i < cell_indices.size(); ++i) {
    adjoint_flux_map[cell_indices[i]] = adjoint_flux_values[i];

    if (adjoint_has_group_data_) {
      const int cell = cell_indices[i];
      std::vector<double> groups(adjoint_n_groups_, 0.0);
      size_t offset = static_cast<size_t>(i) * adjoint_n_groups_;
      for (int g = 0; g < adjoint_n_groups_; ++g) {
        groups[g] = adjoint_group_values[offset + g];
      }
      adjoint_group_map_[cell] = std::move(groups);
    }
  }

  file_close(file_id);
}

//------------------------------------------------------------------------------

void BetaEffective::validate_group_metadata()
{
  const bool flux_multi = flux_has_group_data_ && flux_n_groups_ > 1;
  const bool adjoint_multi = adjoint_has_group_data_ && adjoint_n_groups_ > 1;

  if (!flux_multi || !adjoint_multi) {
    fatal_error(
      "Multi-group β_eff now requires both flux and adjoint HDF5 files to "
      "provide matching group-resolved spectra."
      " Please regenerate the inputs with multi-group tallies enabled.");
  }

  if (flux_n_groups_ != adjoint_n_groups_) {
    fatal_error("Flux and adjoint files provide different n_groups; cannot "
                "perform multi-group β_eff computation.");
  }

  n_energy_groups_ = flux_n_groups_;
  energy_edges_common_.clear();

  if (!flux_energy_edges_.empty() && !adjoint_energy_edges_.empty()) {
    if (flux_energy_edges_.size() != adjoint_energy_edges_.size()) {
      fatal_error("Flux and adjoint energy grids have different sizes."
                  " Please ensure both files share identical edges.");
    }

    const double tol = 1e-8;
    for (size_t i = 0; i < flux_energy_edges_.size(); ++i) {
      double a = flux_energy_edges_[i];
      double b = adjoint_energy_edges_[i];
      if (std::abs(a - b) > tol * std::max(1.0, std::abs(a))) {
        fatal_error("Flux/adjoint energy grid mismatch detected at edge " +
                    std::to_string(i));
      }
    }
    energy_edges_common_ = flux_energy_edges_;
  } else if (!flux_energy_edges_.empty()) {
    energy_edges_common_ = flux_energy_edges_;
  } else if (!adjoint_energy_edges_.empty()) {
    energy_edges_common_ = adjoint_energy_edges_;
  }

  if (flux_group_map_.empty() || adjoint_group_map_.empty()) {
    fatal_error("Flux/adjoint inputs declare multi-group data but contain no "
                "group spectra; ensure tallies were written correctly.");
  }

  std::cout << "  多群 β_eff 已启用 (" << n_energy_groups_ << " 个能群)"
            << std::endl;
}

//------------------------------------------------------------------------------

void BetaEffective::initialize_flux_spectrum_weights()
{
  if (n_energy_groups_ <= 1 || !flux_has_group_data_) {
    flux_collapse_weights_.assign(1, 1.0);
    return;
  }

  flux_collapse_weights_.assign(n_energy_groups_, 0.0);

  for (const auto& [cell_idx, spectrum] : flux_group_map_) {
    if (spectrum.size() != static_cast<size_t>(n_energy_groups_))
      continue;

    for (int g = 0; g < n_energy_groups_; ++g) {
      double value = spectrum[g];
      if (value > 0.0 && std::isfinite(value)) {
        flux_collapse_weights_[g] += value;
      }
    }
  }

  double total = std::accumulate(
    flux_collapse_weights_.begin(), flux_collapse_weights_.end(), 0.0);

  if (total <= 0.0) {
    double uniform = 1.0 / static_cast<double>(n_energy_groups_);
    std::fill(
      flux_collapse_weights_.begin(), flux_collapse_weights_.end(), uniform);
    std::cout << "  通量谱权重不可用，使用均匀折合权重" << std::endl;
  } else {
    for (double& value : flux_collapse_weights_) {
      value /= total;
    }
    std::cout << "  通量加权折合已启用 (" << flux_group_map_.size()
              << " 个网格单元)" << std::endl;
  }
}

//------------------------------------------------------------------------------

int BetaEffective::group_index_from_energy(double energy_eV) const
{
  if (n_energy_groups_ <= 1 || energy_edges_common_.size() < 2) {
    return 0;
  }

  if (energy_eV <= energy_edges_common_.front()) {
    return 0;
  }
  if (energy_eV >= energy_edges_common_.back()) {
    return n_energy_groups_ - 1;
  }

  auto it = std::upper_bound(
    energy_edges_common_.begin(), energy_edges_common_.end(), energy_eV);
  int idx =
    static_cast<int>(std::distance(energy_edges_common_.begin(), it)) - 1;
  if (idx < 0)
    idx = 0;
  if (idx >= n_energy_groups_)
    idx = n_energy_groups_ - 1;
  return idx;
}

//------------------------------------------------------------------------------

size_t BetaEffective::delayed_offset(int energy_group, int delayed_group) const
{
  return static_cast<size_t>(energy_group) * N_DELAYED_GROUPS + delayed_group;
}

//------------------------------------------------------------------------------

double BetaEffective::get_beta_i(int group) const
{
  if (group < 0 || group >= 8) {
    fatal_error("Invalid delayed group index: " + std::to_string(group) +
                " (must be 0-7)");
  }
  return beta_i_[group];
}

//------------------------------------------------------------------------------

void BetaEffective::write_to_file(const std::string& filename) const
{
  hid_t file_id = file_open(filename, 'w');

  // 写入文件元数据属性
  write_attribute(file_id, "filetype", "beta_effective");
  write_attribute(file_id, "version", "1.0");
  write_attribute(
    file_id, "description", "Effective delayed neutron fraction (β_eff)");

  // ========== 主要结果数据集 ==========
  // beta_i [8] - 各组 β_i,eff (Phase 3: 8组)
  std::vector<double> beta_i_vec(beta_i_.begin(), beta_i_.end());
  write_dataset(file_id, "beta_i", beta_i_vec);

  // beta_total - 总 β_eff
  std::vector<double> beta_total_vec = {beta_total_};
  write_dataset(file_id, "beta_total", beta_total_vec);

  // uncertainty [8] - Phase 1暂时为零,预留给Phase 4统计不确定度
  std::vector<double> uncertainty(8, 0.0); // Phase 3: 8组
  write_dataset(file_id, "uncertainty", uncertainty);

  // ========== metadata/ 组 ==========
  hid_t metadata_group =
    H5Gcreate(file_id, "metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  write_attribute(metadata_group, "flux_file", "flux_mesh.h5");
  write_attribute(metadata_group, "adjoint_flux_file", "adjoint_flux.h5");

  // 记录计算时间(秒级时间戳)
  auto now = std::chrono::system_clock::now();
  auto timestamp = std::chrono::system_clock::to_time_t(now);
  write_attribute(
    metadata_group, "computation_time", static_cast<int64_t>(timestamp));

  // 额外的元数据
  std::string nuclear_data_source =
    "Phase 3 - Dynamic extraction from OpenMC library";
  write_attribute(metadata_group, "nuclear_data_source", nuclear_data_source);
  write_attribute(metadata_group, "energy_groups", n_energy_groups_);
  if (!energy_edges_common_.empty()) {
    write_dataset(metadata_group, "energy_edges", energy_edges_common_);
  }
  write_attribute(metadata_group, "beta_eff_mode", "MATERIAL_DEPENDENT");

  // Phase 3: 添加参考能量和温度
  write_attribute(metadata_group, "ref_energy_ev", ref_energy_);
  write_attribute(metadata_group, "ref_temperature_k", ref_temperature_);
  write_attribute(metadata_group, "n_delayed_groups", 8); // Phase 3: 支持8组

  write_dataset(metadata_group, "grid_shape", grid_shape_);
  write_attribute(metadata_group, "grid_pitch", grid_pitch_);

  H5Gclose(metadata_group);

  // ========== diagnostics/ 组 ==========
  hid_t diag_group =
    H5Gcreate(file_id, "diagnostics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // numerator [8] - 各组分子项 (Phase 3: 8组)
  std::vector<double> numerators_vec(numerators_.begin(), numerators_.end());
  write_dataset(diag_group, "numerator", numerators_vec);

  // denominator - 分母项(单个值)
  std::vector<double> denom_vec = {denominator_};
  write_dataset(diag_group, "denominator", denom_vec);

  // normalization_factor - 归一化因子(体积元)
  write_attribute(diag_group, "normalization_factor",
    grid_pitch_ * grid_pitch_ * grid_pitch_);

  // 输出材料信息
  std::vector<int> material_ids;
  std::vector<double> material_nu_total;
  std::vector<double> material_sigma_f;
  std::vector<int> material_cell_counts; // Phase 2.2: 每个材料的单元数

  // 统计材料分布
  std::unordered_map<int, int> mat_counts;
  for (const auto& [cell_idx, mat_id] : cell_to_material_) {
    mat_counts[mat_id]++;
  }

  for (const auto& mat_data : unique_materials_) {
    material_ids.push_back(mat_data.material_id);
    material_nu_total.push_back(mat_data.nu_total);
    material_sigma_f.push_back(mat_data.sigma_f);
    material_cell_counts.push_back(mat_counts[mat_data.material_id]);
  }

  if (!material_ids.empty()) {
    write_dataset(diag_group, "material_ids", material_ids);
    write_dataset(diag_group, "material_nu_total", material_nu_total);
    write_dataset(diag_group, "material_sigma_f", material_sigma_f);
    write_dataset(diag_group, "material_cell_counts", material_cell_counts);
    write_attribute(
      diag_group, "n_unique_materials", static_cast<int>(material_ids.size()));
    write_attribute(diag_group, "total_fissionable_cells",
      static_cast<int>(cell_to_material_.size()));
  }

  H5Gclose(diag_group);
  file_close(file_id);

  std::cout << "  β_eff = " << std::fixed << std::setprecision(5) << beta_total_
            << " (基于" << unique_materials_.size() << "种裂变材料)"
            << std::endl;
}
//==============================================================================
// Material-Dependent Methods
//==============================================================================

//------------------------------------------------------------------------------
// 材料感知的缓发分子（第 i 组）：
//   N_i = Σ_cell ΔV × [Σ_{g'} φ*_{g'} χ_{d,i,g'}] × [Σ_g ν_{d,i,g} Σ_{f,g} φ_g]
// 出生能量 E(g') 和诱发裂变能量 E'(g) 是两个独立积分变量，
// 不能在同一群内逐项相乘。若无多群信息则退化为单群近似。
//------------------------------------------------------------------------------
double BetaEffective::compute_delayed_numerator_material(int group,
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  const bool use_multi_group =
    n_energy_groups_ > 1 && flux_has_group_data_ && adjoint_has_group_data_ &&
    !flux_group_map_.empty() && !adjoint_group_map_.empty();

  // 使用 std::once_flag 保证线程安全的一次性输出
  static std::once_flag numerator_debug_flag;
  std::call_once(numerator_debug_flag, [&]() {
    std::cout << "  [debug] use_multi_group = "
              << (use_multi_group ? "true" : "false") << std::endl;
  });

  std::vector<double> terms;
  terms.reserve(use_multi_group ? flux.size() * n_energy_groups_ : flux.size());

  for (const auto& [cell_idx, phi] : flux) {
    auto it_adj_total = adjoint_flux.find(cell_idx);
    if (it_adj_total == adjoint_flux.end())
      continue;

    auto it_data = cell_nuclear_data_.find(cell_idx);
    if (it_data == cell_nuclear_data_.end() || !it_data->second.is_fissionable)
      continue;

    const auto& nuc_data = it_data->second;

    if (!use_multi_group) {
      double phi_star = it_adj_total->second;
      double term =
        phi_star * nuc_data.nu_delayed[group] * nuc_data.sigma_f * phi * volume;
      terms.push_back(term);
      continue;
    }

    auto it_flux_groups = flux_group_map_.find(cell_idx);
    auto it_adj_groups = adjoint_group_map_.find(cell_idx);
    if (it_flux_groups == flux_group_map_.end() ||
        it_adj_groups == adjoint_group_map_.end())
      continue;

    const auto& flux_groups = it_flux_groups->second;
    const auto& adjoint_groups = it_adj_groups->second;
    if (flux_groups.size() != static_cast<size_t>(n_energy_groups_) ||
        adjoint_groups.size() != static_cast<size_t>(n_energy_groups_)) {
      continue;
    }

    if (nuc_data.sigma_f_groups.size() !=
          static_cast<size_t>(n_energy_groups_) ||
        nuc_data.nu_delayed_groups.size() !=
          static_cast<size_t>(n_energy_groups_) * N_DELAYED_GROUPS ||
        nuc_data.chi_delayed_groups.size() !=
          static_cast<size_t>(n_energy_groups_) * N_DELAYED_GROUPS) {
      fatal_error("Material " + nuc_data.material_name +
                  " missing multi-group delayed data.");
    }

    // 出生能量积分: A_{d,i} = Σ_{g'} φ*_{g'} × χ_{d,i,g'}
    double adj_chi_delayed = 0.0;
    for (int gp = 0; gp < n_energy_groups_; ++gp) {
      adj_chi_delayed += adjoint_groups[gp] *
                         nuc_data.chi_delayed_groups[delayed_offset(gp, group)];
    }

    // 诱发裂变能量积分: F_{d,i} = Σ_g ν_{d,i,g} × Σ_{f,g} × φ_g
    double delayed_fission_source = 0.0;
    for (int g = 0; g < n_energy_groups_; ++g) {
      delayed_fission_source +=
        nuc_data.nu_delayed_groups[delayed_offset(g, group)] *
        nuc_data.sigma_f_groups[g] * flux_groups[g];
    }

    // 单元贡献 = ΔV × A_{d,i} × F_{d,i}
    double term = volume * adj_chi_delayed * delayed_fission_source;
    if (term != 0.0) {
      terms.push_back(term);
    }
  }

  return kahan_sum(terms);
}

//------------------------------------------------------------------------------

// 分母（瞬发项）：
//   D = Σ_cell ΔV × [Σ_{g'} φ*_{g'} χ_p,{g'}] × [Σ_g ν_g Σ_{f,g} φ_g]
// 出生能量 E(g') 和诱发裂变能量 E'(g) 是两个独立积分变量，
// 不能在同一群内逐项相乘。
double BetaEffective::compute_denominator_material(
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  std::vector<double> terms;
  terms.reserve(flux.size() * n_energy_groups_);

  for (const auto& [cell_idx, phi] : flux) {
    auto it_adj_total = adjoint_flux.find(cell_idx);
    if (it_adj_total == adjoint_flux.end())
      continue;

    auto it_data = cell_nuclear_data_.find(cell_idx);
    if (it_data == cell_nuclear_data_.end() || !it_data->second.is_fissionable)
      continue;

    const auto& nuc_data = it_data->second;

    auto it_flux_groups = flux_group_map_.find(cell_idx);
    auto it_adj_groups = adjoint_group_map_.find(cell_idx);
    if (it_flux_groups == flux_group_map_.end() ||
        it_adj_groups == adjoint_group_map_.end())
      continue;

    const auto& flux_groups = it_flux_groups->second;
    const auto& adjoint_groups = it_adj_groups->second;
    if (flux_groups.size() != static_cast<size_t>(n_energy_groups_) ||
        adjoint_groups.size() != static_cast<size_t>(n_energy_groups_)) {
      fatal_error("单元 " + std::to_string(cell_idx) + " 缺少多群通量数据。");
    }

    if (nuc_data.sigma_f_groups.size() !=
          static_cast<size_t>(n_energy_groups_) ||
        nuc_data.nu_total_groups.size() !=
          static_cast<size_t>(n_energy_groups_) ||
        nuc_data.chi_prompt_groups.size() !=
          static_cast<size_t>(n_energy_groups_)) {
      fatal_error("材料 " + nuc_data.material_name + " 缺少多群核数据。");
    }

    // 出生能量积分: A_p = Σ_{g'} φ*_{g'} × χ_{p,g'}
    double adj_chi_prompt = 0.0;
    for (int gp = 0; gp < n_energy_groups_; ++gp) {
      adj_chi_prompt += adjoint_groups[gp] * nuc_data.chi_prompt_groups[gp];
    }

    // 诱发裂变能量积分: F = Σ_g ν_g × Σ_{f,g} × φ_g
    double fission_source = 0.0;
    for (int g = 0; g < n_energy_groups_; ++g) {
      fission_source += nuc_data.nu_total_groups[g] *
                        nuc_data.sigma_f_groups[g] * flux_groups[g];
    }

    // 单元贡献 = ΔV × A_p × F
    double term = volume * adj_chi_prompt * fission_source;
    if (term != 0.0) {
      terms.push_back(term);
    }
  }

  return kahan_sum(terms);
}

//------------------------------------------------------------------------------
// 计算不含 chi 的分母（用于诊断对比）
// 公式: D = Σ_cell Σ_g φ*_g × ν_total_g × Σ_f,g × φ_g × ΔV
// 这是另一种常见的 β_eff 分母定义
//------------------------------------------------------------------------------
double BetaEffective::compute_denominator_without_chi(
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  std::vector<double> terms;
  terms.reserve(flux.size() * n_energy_groups_);

  for (const auto& [cell_idx, phi] : flux) {
    auto it_adj_total = adjoint_flux.find(cell_idx);
    if (it_adj_total == adjoint_flux.end())
      continue;

    auto it_data = cell_nuclear_data_.find(cell_idx);
    if (it_data == cell_nuclear_data_.end() || !it_data->second.is_fissionable)
      continue;

    const auto& nuc_data = it_data->second;

    auto it_flux_groups = flux_group_map_.find(cell_idx);
    auto it_adj_groups = adjoint_group_map_.find(cell_idx);
    if (it_flux_groups == flux_group_map_.end() ||
        it_adj_groups == adjoint_group_map_.end())
      continue;

    const auto& flux_groups = it_flux_groups->second;
    const auto& adjoint_groups = it_adj_groups->second;
    if (flux_groups.size() != static_cast<size_t>(n_energy_groups_) ||
        adjoint_groups.size() != static_cast<size_t>(n_energy_groups_))
      continue;

    if (nuc_data.sigma_f_groups.size() !=
          static_cast<size_t>(n_energy_groups_) ||
        nuc_data.nu_total_groups.size() !=
          static_cast<size_t>(n_energy_groups_))
      continue;

    for (int g = 0; g < n_energy_groups_; ++g) {
      double phi_g = flux_groups[g];
      double phi_star_g = adjoint_groups[g];
      if (phi_g == 0.0 || phi_star_g == 0.0)
        continue;

      double sigma_f_g = nuc_data.sigma_f_groups[g];
      double nu_total_g = nuc_data.nu_total_groups[g];

      // 不含 chi 的分母项
      double term = phi_star_g * nu_total_g * sigma_f_g * phi_g * volume;
      terms.push_back(term);
    }
  }

  return kahan_sum(terms);
}

//------------------------------------------------------------------------------
// 从 OpenMC 核数据库动态提取核参数
//------------------------------------------------------------------------------
MaterialNuclearData BetaEffective::extract_material_nuclear_data(
  int material_id) const
{
  // 校验输入参数
  const int n_groups = n_energy_groups_;
  if (n_groups <= 1) {
    fatal_error("材料数据提取需要多群输入。");
  }

  if (energy_edges_common_.size() != static_cast<size_t>(n_groups + 1)) {
    fatal_error("多群 β_eff 需要 n_groups + 1 个能量边界。");
  }

  // 获取通量折合权重
  std::vector<double> collapse_weights(
    n_groups, 1.0 / static_cast<double>(n_groups));
  if (flux_collapse_weights_.size() == static_cast<size_t>(n_groups)) {
    collapse_weights = flux_collapse_weights_;
  }

  // 使用 MaterialNuclearDataExtractor 提取核数据
  MaterialNuclearDataExtractor extractor(n_groups, energy_edges_common_,
    collapse_weights, ref_energy_, ref_temperature_);

  static std::once_flag extract_log_flag;
  std::call_once(extract_log_flag,
    [&]() { std::cout << "  使用材料核数据提取器" << std::endl; });

  return extractor.extract(material_id);
}

//------------------------------------------------------------------------------

void BetaEffective::build_cell_material_map(
  const std::unordered_map<int, double>& flux)
{
  std::cout << "  构建单元-材料映射 (" << flux.size() << " 个网格单元)..."
            << std::endl;

  cell_to_material_.clear();
  cell_nuclear_data_.clear();
  n_heterogeneous_cells_ = 0;

  // 检查材料库是否可用
  if (model::materials.empty()) {
    warning("Material library is empty! Cannot perform geometry query.");
    return;
  }

  std::cout << "  材料库: " << model::materials.size() << " 种材料"
            << std::endl;

  // ===== 在并行区域之前预先提取所有裂变材料的核数据 =====
  // 这是线程安全的关键：所有核数据库访问都在单线程中完成
  std::unordered_map<int, MaterialNuclearData> global_material_cache;
  std::cout << "  预提取裂变材料核数据..." << std::endl;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->fissionable()) {
      int mat_id = mat_ptr->id();
      global_material_cache[mat_id] = extract_material_nuclear_data(mat_id);
    }
  }
  std::cout << "  已提取 " << global_material_cache.size() << " 种裂变材料数据"
            << std::endl;

  // 线程局部缓存和统计
  struct ThreadLocalData {
    std::unordered_map<int, MaterialNuclearData> cell_nuclear_data;
    std::unordered_map<int, int> cell_to_material;
    std::unordered_map<int, int> material_hit_counts;
    int fissionable_count = 0;
    int non_fissionable_count = 0;
    int void_count = 0;
    int geometry_failed_count = 0;
    int heterogeneous_count = 0;
  };

  std::vector<ThreadLocalData> thread_data;

  // 将 flux map 转换为 vector 以便并行化 (在并行区域外部创建)
  std::vector<std::pair<int, double>> flux_vec(flux.begin(), flux.end());

  // 预先确定线程数并分配存储
  int num_threads_max = 1;
#ifdef _OPENMP
  num_threads_max = omp_get_max_threads();
#endif
  thread_data.resize(num_threads_max);

#pragma omp parallel
  {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif

    auto& local_data = thread_data[thread_id];
    GeometryState geom; // 线程私有几何状态

#pragma omp for schedule(dynamic, 100)
    for (size_t idx = 0; idx < flux_vec.size(); ++idx) {
      int cell_idx = flux_vec[idx].first;
      double flux_val = flux_vec[idx].second;

      // 计算网格单元边界
      Position center_pos = grid_index_to_position(cell_idx);
      Position cell_lower = center_pos - Position {grid_pitch_ / 2,
                                           grid_pitch_ / 2, grid_pitch_ / 2};
      Position cell_upper = center_pos + Position {grid_pitch_ / 2,
                                           grid_pitch_ / 2, grid_pitch_ / 2};

      // 多点采样
      std::unordered_map<int, int> material_counts;
      int fissionable_hits = 0;     // 裂变材料命中次数
      int non_fissionable_hits = 0; // 非裂变材料命中次数
      int void_hits = 0;            // 空白区域命中次数

      // 对每个网格单元执行 1/8/27 点采样来估计材料体积分数。
      // 采样点越多，异质结构估计越精准，但也需要更多 `exhaustive_find_cell`
      // 几何查询。
      for (int i = 0; i < n_sample_points_; ++i) {
        Position sample_pos;

        if (n_sample_points_ == 1) {
          // Phase 2.2: 单点中心采样
          sample_pos = center_pos;
        } else {
          // Phase 2.3: 多点规则网格采样 (8 或 27 点)
          sample_pos =
            sample_cell_point(cell_lower, cell_upper, i, n_sample_points_);
        }

        // 几何查询
        geom.r() = sample_pos;
        geom.u() = Direction {0.0, 0.0, 1.0};

        if (!exhaustive_find_cell(geom)) {
          continue; // 几何查询失败
        }

        // 获取材料索引
        int cell_index = geom.lowest_coord().cell();
        Cell* cell = model::cells[cell_index].get();
        int instance = geom.cell_instance();
        int material_idx = cell->material(instance);

        if (material_idx == MATERIAL_VOID) {
          void_hits++; // 统计空白区域
          continue;
        }

        // 验证材料索引
        if (material_idx < 0 ||
            material_idx >= static_cast<int>(model::materials.size())) {
          continue;
        }

        const auto& material = model::materials[material_idx];
        if (!material) {
          void_hits++; // 无效材料视为空白
          continue;
        }

        if (!material->fissionable()) {
          non_fissionable_hits++; // 统计非裂变材料(慢化剂、结构等)
          continue;
        }

        // 裂变材料: 参与体积分数计算
        int material_id = material->id();
        material_counts[material_id]++;
        fissionable_hits++;
      }

      // 处理采样结果
      if (fissionable_hits == 0) {
        // 没有裂变材料命中
        if (void_hits > 0 && non_fissionable_hits == 0) {
          // 纯空白单元
          local_data.void_count++;
        } else {
          // 仅含非裂变材料的单元
          local_data.non_fissionable_count++;
        }
        continue;
      }

      // 检查是否异质(多种裂变材料)
      bool is_heterogeneous = (material_counts.size() > 1);
      if (is_heterogeneous) {
        local_data.heterogeneous_count++;
      }

      MaterialNuclearData cell_data;

      if (material_counts.size() == 1) {
        // 同质单元: 直接从全局缓存获取材料数据（线程安全，只读）
        int mat_id = material_counts.begin()->first;
        auto it = global_material_cache.find(mat_id);
        if (it != global_material_cache.end()) {
          cell_data = it->second;
        }
        local_data.cell_to_material[cell_idx] = mat_id;
      } else {
        // 异质单元: 使用全局缓存计算加权核数据（线程安全，只读）
        cell_data = compute_weighted_nuclear_data(
          material_counts, fissionable_hits, global_material_cache);
        local_data.cell_to_material[cell_idx] = -1; // 标记为多材料
      }

      local_data.cell_nuclear_data[cell_idx] = cell_data;
      local_data.fissionable_count++;

      // 统计材料分布
      for (const auto& [mat_id, count] : material_counts) {
        local_data.material_hit_counts[mat_id] += count;
      }
    }
  } // end parallel region

  // 合并线程结果
  std::unordered_map<int, int> material_hit_counts;
  int fissionable_count = 0;
  int non_fissionable_count = 0;
  int void_count = 0;
  int geometry_failed_count = 0;

  for (const auto& local : thread_data) {
    // 合并单元数据
    for (const auto& [cell_idx, data] : local.cell_nuclear_data) {
      cell_nuclear_data_[cell_idx] = data;
    }
    for (const auto& [cell_idx, mat_id] : local.cell_to_material) {
      cell_to_material_[cell_idx] = mat_id;
    }

    // 合并统计
    for (const auto& [mat_id, count] : local.material_hit_counts) {
      material_hit_counts[mat_id] += count;
    }
    fissionable_count += local.fissionable_count;
    non_fissionable_count += local.non_fissionable_count;
    void_count += local.void_count;
    geometry_failed_count += local.geometry_failed_count;
    n_heterogeneous_cells_ += local.heterogeneous_count;
  }

  // 构建唯一材料列表（使用全局缓存）
  unique_materials_.clear();
  for (const auto& [mat_id, data] : global_material_cache) {
    unique_materials_.push_back(data);
  }

  // 输出统计信息
  std::cout << "  几何查询结果: 裂变单元 " << fissionable_count << ", 非裂变 "
            << non_fissionable_count << ", 空白 " << void_count << ", 材料种类 "
            << unique_materials_.size() << std::endl;
}

//------------------------------------------------------------------------------

std::array<double, 3> BetaEffective::mesh_index_to_position(int mesh_idx) const
{
  // 将一维索引转换为三维坐标
  int nx = grid_shape_[0];
  int ny = grid_shape_[1];
  int nz = grid_shape_[2];

  int iz = mesh_idx / (nx * ny);
  int remainder = mesh_idx % (nx * ny);
  int iy = remainder / nx;
  int ix = remainder % nx;

  // 计算单元中心位置
  // 假设网格从 lower_left 开始，每个单元大小为 grid_pitch_
  double x = grid_lower_left_[0] + (ix + 0.5) * grid_pitch_;
  double y = grid_lower_left_[1] + (iy + 0.5) * grid_pitch_;
  double z = grid_lower_left_[2] + (iz + 0.5) * grid_pitch_;

  return {x, y, z};
}

//------------------------------------------------------------------------------

std::array<int, 3> BetaEffective::get_grid_indices(int cell_idx) const
{
  int nx = grid_shape_[0];
  int ny = grid_shape_[1];
  int nz = grid_shape_[2];

  // 从一维索引反推三维索引 (row-major order: z,y,x)
  int iz = cell_idx / (nx * ny);
  int remainder = cell_idx % (nx * ny);
  int iy = remainder / nx;
  int ix = remainder % nx;

  return {ix, iy, iz};
}

//------------------------------------------------------------------------------

Position BetaEffective::grid_index_to_position(int cell_idx) const
{
  // 获取三维网格索引
  auto [ix, iy, iz] = get_grid_indices(cell_idx);

  // 计算网格单元中心的全局坐标
  // 单元中心 = lower_left + (index + 0.5) * pitch
  double x = grid_lower_left_[0] + (ix + 0.5) * grid_pitch_;
  double y = grid_lower_left_[1] + (iy + 0.5) * grid_pitch_;
  double z = grid_lower_left_[2] + (iz + 0.5) * grid_pitch_;

  return Position {x, y, z};
}

//------------------------------------------------------------------------------

Position BetaEffective::sample_cell_point(const Position& lower,
  const Position& upper, int point_idx, int n_points) const
{
  if (n_points == 8) {
    // 8 个立方体顶点 (000, 100, 010, 110, 001, 101, 011, 111)
    const double offsets[8][3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
      {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

    if (point_idx < 0 || point_idx >= 8) {
      fatal_error("Invalid corner index: " + std::to_string(point_idx) +
                  " (must be 0-7)");
    }

    return Position {lower.x + offsets[point_idx][0] * (upper.x - lower.x),
      lower.y + offsets[point_idx][1] * (upper.y - lower.y),
      lower.z + offsets[point_idx][2] * (upper.z - lower.z)};

  } else if (n_points == 27) {
    // 3×3×3 规则网格采样
    if (point_idx < 0 || point_idx >= 27) {
      fatal_error("Invalid point index: " + std::to_string(point_idx) +
                  " (must be 0-26 for 27-point sampling)");
    }

    // 从点索引计算 3D 网格坐标 (i,j,k)
    int iz = point_idx / 9; // 0, 1, 2
    int remainder = point_idx % 9;
    int iy = remainder / 3; // 0, 1, 2
    int ix = remainder % 3; // 0, 1, 2

    // 计算相对位置 (0, 0.5, 1)
    double fx = ix * 0.5;
    double fy = iy * 0.5;
    double fz = iz * 0.5;

    return Position {lower.x + fx * (upper.x - lower.x),
      lower.y + fy * (upper.y - lower.y), lower.z + fz * (upper.z - lower.z)};

  } else {
    fatal_error("Unsupported number of sample points: " +
                std::to_string(n_points) + " (must be 8 or 27)");
    return Position {0, 0, 0}; // 永远不会到达
  }
}

//------------------------------------------------------------------------------

MaterialNuclearData BetaEffective::compute_weighted_nuclear_data(
  const std::unordered_map<int, int>& material_counts, int total_samples,
  const std::unordered_map<int, MaterialNuclearData>& material_cache) const
{
  // 将异质单元折算为单一“有效”材料。调用者提供每个材料的命中次数，函数
  // 用体积分数或裂变反应率作为权重对 Σ_f、ν、χ 等量做混合，使后续公式仍可
  // 使用单一 MaterialNuclearData 接口而无需特殊分支。
  MaterialNuclearData weighted_data;
  weighted_data.material_id = -1; // 标记为多材料单元
  weighted_data.material_name = "mixed";
  weighted_data.is_fissionable = true;

  const int n_groups = std::max(1, n_energy_groups_);
  const size_t delayed_size = static_cast<size_t>(n_groups) * N_DELAYED_GROUPS;

  std::vector<double> accum_sigma_groups(n_groups, 0.0);
  std::vector<double> accum_nu_total_groups(n_groups, 0.0);
  std::vector<double> accum_nu_prompt_groups(n_groups, 0.0);
  std::vector<double> accum_nu_delayed_groups(delayed_size, 0.0);
  std::vector<double> accum_chi_prompt(n_groups, 0.0);
  std::vector<double> accum_chi_delayed(delayed_size, 0.0);

  auto group_value = [&](const std::vector<double>& values, double fallback,
                       int g) {
    if (values.size() == static_cast<size_t>(n_groups)) {
      return values[g];
    }
    return fallback;
  };

  auto delayed_value = [&](const std::vector<double>& values,
                         const std::array<double, 8>& fallback, int g, int d) {
    if (values.size() == delayed_size) {
      return values[delayed_offset(g, d)];
    }
    return fallback[d];
  };

  auto accumulate_group_data = [&](const MaterialNuclearData& mat_data,
                                 double weight) {
    for (int g = 0; g < n_groups; ++g) {
      accum_sigma_groups[g] +=
        weight * group_value(mat_data.sigma_f_groups, mat_data.sigma_f, g);
      accum_nu_total_groups[g] +=
        weight * group_value(mat_data.nu_total_groups, mat_data.nu_total, g);
      accum_nu_prompt_groups[g] +=
        weight * group_value(mat_data.nu_prompt_groups, mat_data.nu_prompt, g);

      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        accum_nu_delayed_groups[delayed_offset(g, d)] +=
          weight *
          delayed_value(mat_data.nu_delayed_groups, mat_data.nu_delayed, g, d);
        double chi_delayed_val = delayed_value(
          mat_data.chi_delayed_groups, mat_data.chi_delayed, g, d);
        accum_chi_delayed[delayed_offset(g, d)] += weight * chi_delayed_val;
      }

      double chi_prompt_val =
        group_value(mat_data.chi_prompt_groups, mat_data.chi_prompt, g);
      accum_chi_prompt[g] += weight * chi_prompt_val;
    }
  };

  auto normalize_distribution = [](std::vector<double>& dist) {
    if (dist.empty())
      return;
    double total = std::accumulate(dist.begin(), dist.end(), 0.0);
    if (total <= 0.0) {
      double uniform = 1.0 / static_cast<double>(dist.size());
      std::fill(dist.begin(), dist.end(), uniform);
      return;
    }
    for (double& val : dist) {
      val /= total;
    }
  };

  if (weighting_mode_ == WeightingMode::VOLUME_WEIGHTED) {
    // ===== 体积分数加权 =====
    for (const auto& [mat_id, count] : material_counts) {
      double fraction = static_cast<double>(count) / total_samples;

      // 从缓存获取材料核数据（线程安全）
      auto it = material_cache.find(mat_id);
      if (it == material_cache.end())
        continue;
      const MaterialNuclearData& mat_data = it->second;

      // 简单加权平均
      weighted_data.nu_total += fraction * mat_data.nu_total;
      weighted_data.nu_prompt += fraction * mat_data.nu_prompt;
      weighted_data.sigma_f += fraction * mat_data.sigma_f;

      for (int g = 0; g < 8; ++g) { // Phase 3: 8组
        weighted_data.nu_delayed[g] += fraction * mat_data.nu_delayed[g];
      }

      accumulate_group_data(mat_data, fraction);
    }

  } else {
    // ===== 反应率加权 (物理更准确) =====
    double total_fission_density = 0.0;
    double weighted_nu_total = 0.0;
    double weighted_nu_prompt = 0.0;
    std::array<double, 8> weighted_nu_delayed = {}; // Phase 3: 8组

    for (const auto& [mat_id, count] : material_counts) {
      double fraction = static_cast<double>(count) / total_samples;

      // 从缓存获取材料核数据（线程安全）
      auto it = material_cache.find(mat_id);
      if (it == material_cache.end())
        continue;
      const MaterialNuclearData& mat_data = it->second;

      // 按裂变反应率加权: fraction × Σ_f
      double fission_contribution = fraction * mat_data.sigma_f;

      total_fission_density += fission_contribution;
      weighted_data.sigma_f += fission_contribution;

      // 累积加权求和
      weighted_nu_total += fission_contribution * mat_data.nu_total;
      weighted_nu_prompt += fission_contribution * mat_data.nu_prompt;

      for (int g = 0; g < 8; ++g) { // Phase 3: 8组
        weighted_nu_delayed[g] += fission_contribution * mat_data.nu_delayed[g];
      }

      accumulate_group_data(mat_data, fission_contribution);
    }

    // 归一化：除以总裂变密度
    if (total_fission_density > 0.0) {
      weighted_data.nu_total = weighted_nu_total / total_fission_density;
      weighted_data.nu_prompt = weighted_nu_prompt / total_fission_density;
      for (int g = 0; g < 8; ++g) { // Phase 3: 8组
        weighted_data.nu_delayed[g] =
          weighted_nu_delayed[g] / total_fission_density;
      }
      // sigma_f 已经是加权后的宏观截面，不需要归一化

      for (double& value : accum_sigma_groups) {
        value /= total_fission_density;
      }
      for (double& value : accum_nu_total_groups) {
        value /= total_fission_density;
      }
      for (double& value : accum_nu_prompt_groups) {
        value /= total_fission_density;
      }
      for (double& value : accum_nu_delayed_groups) {
        value /= total_fission_density;
      }
      for (double& value : accum_chi_prompt) {
        value /= total_fission_density;
      }
      for (double& value : accum_chi_delayed) {
        value /= total_fission_density;
      }
    }
  }

  // 规范化 χ 分布
  normalize_distribution(accum_chi_prompt);
  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    std::vector<double> delayed_slice(n_groups, 0.0);
    for (int g = 0; g < n_groups; ++g) {
      delayed_slice[g] = accum_chi_delayed[delayed_offset(g, d)];
    }
    normalize_distribution(delayed_slice);
    for (int g = 0; g < n_groups; ++g) {
      accum_chi_delayed[delayed_offset(g, d)] = delayed_slice[g];
    }
  }

  weighted_data.sigma_f_groups = accum_sigma_groups;
  weighted_data.nu_total_groups = accum_nu_total_groups;
  weighted_data.nu_prompt_groups = accum_nu_prompt_groups;
  weighted_data.nu_delayed_groups = accum_nu_delayed_groups;
  weighted_data.chi_prompt_groups = accum_chi_prompt;
  weighted_data.chi_delayed_groups = accum_chi_delayed;

  auto ensure_group_values = [&](std::vector<double>& target, double fallback) {
    if (target.empty())
      return;
    bool all_zero = std::all_of(target.begin(), target.end(),
      [](double v) { return std::abs(v) < 1e-16; });
    if (all_zero) {
      std::fill(target.begin(), target.end(), fallback);
    }
  };

  ensure_group_values(weighted_data.sigma_f_groups, weighted_data.sigma_f);
  ensure_group_values(weighted_data.nu_total_groups, weighted_data.nu_total);
  ensure_group_values(weighted_data.nu_prompt_groups, weighted_data.nu_prompt);

  if (!weighted_data.nu_delayed_groups.empty()) {
    bool delayed_zero = std::all_of(weighted_data.nu_delayed_groups.begin(),
      weighted_data.nu_delayed_groups.end(),
      [](double v) { return std::abs(v) < 1e-16; });
    if (delayed_zero) {
      for (int g = 0; g < n_groups; ++g) {
        for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
          weighted_data.nu_delayed_groups[delayed_offset(g, d)] =
            weighted_data.nu_delayed[d];
        }
      }
    }
  }

  return weighted_data;
}

//------------------------------------------------------------------------------

double BetaEffective::kahan_sum(const std::vector<double>& values) const
{
  if (values.empty())
    return 0.0;

  // Kahan 求和算法（补偿求和）
  double sum = 0.0;
  double compensation = 0.0; // 累积的舍入误差补偿

  for (double value : values) {
    double y = value - compensation; // 减去上次的误差
    double t = sum + y;              // 加到当前和
    compensation = (t - sum) - y;    // 计算新的舍入误差
    sum = t;                         // 更新和
  }

  return sum;
}

} // namespace openmc