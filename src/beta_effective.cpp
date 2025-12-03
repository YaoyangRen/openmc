#include "openmc/beta_effective.h"

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
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
#include <numeric>
#include <vector>

namespace openmc {

//------------------------------------------------------------------------------
// 核心驱动函数：根据网格通量/共轭通量文件完成 β_eff 计算。
// 算法遵循扰动理论定义：
//   β_i,eff = ∫ φ*(r,E) χ_d,i(E) ν_d,i(E,r) Σ_f(r,E) φ(r,E) dV dE / 分母
// 分母使用瞬发量。实现中分别保留 Phase 1（均匀 U-235 参考核数据）与 Phase 2
// （几何关联、材料相关核数据）的分支，但共享数据读取、网格一致性校验以及
// 输出流程。
//------------------------------------------------------------------------------
void BetaEffective::compute_from_files(const std::string& flux_file,
  const std::string& adjoint_flux_file, const std::string& output_file)
{
  std::cout << std::string(70, '=') << std::endl;

  // 1. 读取正向通量
  std::cout << "[1/4] Reading forward flux from: " << flux_file << std::endl;
  if (!file_exists(flux_file)) {
    fatal_error("Flux file not found: " + flux_file +
                "\nPlease run the forward calculation first.");
  }

  std::unordered_map<int, double> flux;
  std::array<int, 3> flux_shape {};
  double flux_pitch = 0.0;
  read_flux_data(flux_file, flux, flux_shape, flux_pitch);

  std::cout << "  Grid: " << flux_shape[0] << " x " << flux_shape[1] << " x "
            << flux_shape[2] << std::endl;
  std::cout << "  Pitch: " << flux_pitch << " cm" << std::endl;
  std::cout << "  Non-zero cells: " << flux.size() << std::endl;

  // 2. 读取共轭通量 (φ*)，后续作为加权函数将所有单元联系在一起
  std::cout << "\n[2/4] Reading adjoint flux from: " << adjoint_flux_file
            << std::endl;
  if (!file_exists(adjoint_flux_file)) {
    fatal_error("Adjoint flux file not found: " + adjoint_flux_file +
                "\nPlease ensure adjoint flux has been computed.");
  }

  std::unordered_map<int, double> adjoint_flux;
  std::array<int, 3> adjoint_shape {};
  double adjoint_pitch = 0.0;
  read_adjoint_flux_data(
    adjoint_flux_file, adjoint_flux, adjoint_shape, adjoint_pitch);

  std::cout << "  Grid: " << adjoint_shape[0] << " x " << adjoint_shape[1]
            << " x " << adjoint_shape[2] << std::endl;
  std::cout << "  Pitch: " << adjoint_pitch << " cm" << std::endl;
  std::cout << "  Non-zero cells: " << adjoint_flux.size() << std::endl;

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

  std::fill(beta_i_.begin(), beta_i_.end(), 0.0);
  std::fill(numerators_.begin(), numerators_.end(), 0.0);
  denominator_ = 0.0;

  // 计算网格单元体积
  double volume = flux_pitch * flux_pitch * flux_pitch;

  // 根据模式选择计算方法：Phase 1 直接使用常数核数据；Phase 2 需要几何查询
  // 与材料抽样以提取 Σ_f、ν、χ 等空间分布。
  if (mode_ == BetaEffMode::FIXED_U235) {
    // ===== Phase 1: 固定 U-235 核数据 =====
    std::cout << "\n[3/4] Computing β_eff using FIXED U-235 nuclear data"
              << std::endl;
    std::cout << "  Mode: Phase 1 (FIXED_U235)" << std::endl;
    std::cout << "  ν_total = " << NU_TOTAL << std::endl;
    std::cout << "  ν_prompt = " << NU_PROMPT << std::endl;
    std::cout << "  Cell volume = " << volume << " cm³" << std::endl;

    // 计算分母：∫ φ* χ_p ν Σ_f φ dV。Phase 1 不考虑能量依赖，直接使用常数。
    denominator_ = compute_denominator(flux, adjoint_flux, volume);

    if (denominator_ <= 0.0) {
      fatal_error("Denominator is zero or negative! Cannot compute β_eff.");
    }

    std::cout << "  Denominator = " << std::scientific << std::setprecision(6)
              << denominator_ << std::endl;

    // 对每个缓发群计算分子和 β_i,eff (Phase 3: 8组)。在单群近似下仅需
    // 遍历所有非零网格单元。
    std::cout << "\n  Delayed group contributions:" << std::endl;
    std::cout << "  Group    ν_d,i      Numerator        β_i,eff" << std::endl;
    std::cout << "  " << std::string(55, '-') << std::endl;

    for (int i = 0; i < 8; ++i) {
      numerators_[i] = compute_delayed_numerator(i, flux, adjoint_flux, volume);
      beta_i_[i] = numerators_[i] / denominator_;

      std::cout << "    " << (i + 1) << "    " << std::scientific
                << std::setprecision(6) << NU_DELAYED[i] << "   "
                << numerators_[i] << "   " << beta_i_[i] << std::endl;
    }

  } else {
    // ===== Phase 2: 材料相关核数据 =====
    std::cout << "\n[3/4] Computing β_eff using MATERIAL-DEPENDENT nuclear data"
              << std::endl;
    if (n_sample_points_ == 1) {
      std::cout << "  Mode: Phase 2.2 (Single-point geometry query)"
                << std::endl;
    } else {
      std::cout << "  Mode: Phase 2.3 (Multi-point sampling)" << std::endl;
      std::cout << "  Sample points per cell: " << n_sample_points_
                << std::endl;
      std::cout << "  Weighting mode: "
                << (weighting_mode_ == WeightingMode::VOLUME_WEIGHTED
                       ? "Volume-weighted"
                       : "Reaction-rate weighted")
                << std::endl;
    }
    std::cout << "  Cell volume = " << volume << " cm³" << std::endl;

    // 构建单元-材料映射并提取核数据：通过几何查询找出每个网格实际包含的
    // 材料组合，从而得到空间依赖的 Σ_f、ν_t、ν_d、χ。
    build_cell_material_map(flux);

    // 计算分母
    denominator_ = compute_denominator_material(flux, adjoint_flux, volume);

    if (denominator_ <= 0.0) {
      fatal_error("Denominator is zero or negative! Cannot compute β_eff.");
    }

    std::cout << "  Denominator = " << std::scientific << std::setprecision(6)
              << denominator_ << std::endl;

    // 对每个缓发群计算分子和 β_i,eff。若存在多群数据则逐群折合；否则退化
    // 为单群形式。
    std::cout << "\n  Delayed group contributions:" << std::endl;
    std::cout << "  Group    Numerator        β_i,eff" << std::endl;
    std::cout << "  " << std::string(50, '-') << std::endl;

    for (int i = 0; i < 8; ++i) {
      numerators_[i] =
        compute_delayed_numerator_material(i, flux, adjoint_flux, volume);
      beta_i_[i] = numerators_[i] / denominator_;

      std::cout << "    " << (i + 1) << "    " << std::scientific
                << std::setprecision(6) << numerators_[i] << "   " << beta_i_[i]
                << std::endl;
    }
  }

  // 6. 计算总 β_eff
  beta_total_ = std::accumulate(beta_i_.begin(), beta_i_.end(), 0.0);

  std::cout << "  " << std::string(55, '-') << std::endl;
  std::cout << "  Total β_eff = " << std::fixed << std::setprecision(15)
            << beta_total_ << std::endl;

  // 7. 输出结果到文件
  std::cout << "\n[4/4] Writing results to: " << output_file << std::endl;
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

  // 读取网格左下角坐标 (Phase 2 需要)
  if (mode_ == BetaEffMode::MATERIAL_DEPENDENT &&
      object_exists(file_id, "grid_lower_left")) {
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
  read_dataset(file_id, "shape", grid_shape);
  shape = grid_shape;

  read_attribute(file_id, "pitch", pitch);

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
  read_dataset(file_id, "adjoint_flux_values", adjoint_flux_values);

  if (adjoint_n_groups_ > 1 &&
      object_exists(file_id, "adjoint_flux_group_values")) {
    read_dataset(file_id, "adjoint_flux_group_values", adjoint_group_values);
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

  if (flux_multi && adjoint_multi) {
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

    std::cout << "  Multi-group β_eff enabled (" << n_energy_groups_
              << " energy groups)." << std::endl;
  } else {
    if (flux_multi != adjoint_multi) {
      warning("Only one of flux/adjoint inputs contains group data; falling "
              "back to single-group β_eff evaluation.");
    }

    n_energy_groups_ = 1;
    energy_edges_common_.clear();
    flux_group_map_.clear();
    adjoint_group_map_.clear();
    flux_has_group_data_ = false;
    adjoint_has_group_data_ = false;
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

double BetaEffective::compute_delayed_numerator(int group,
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  double sum = 0.0;

  // 遍历所有有正向通量的单元
  for (const auto& [cell_idx, phi] : flux) {
    // 检查该单元是否也有共轭通量
    auto it_adj = adjoint_flux.find(cell_idx);
    if (it_adj == adjoint_flux.end())
      continue;

    double phi_star = it_adj->second;

    // 单能量公式: φ*(r) × χ_d,i × ν_d,i × σ_f × φ(r) × V
    // 使用固定核数据
    sum += phi_star * CHI_DELAYED[group] * NU_DELAYED[group] * SIGMA_F * phi *
           volume;
  }

  return sum;
}

//------------------------------------------------------------------------------

double BetaEffective::compute_denominator(
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  double sum = 0.0;

  // 遍历所有有正向通量的单元
  for (const auto& [cell_idx, phi] : flux) {
    // 检查该单元是否也有共轭通量
    auto it_adj = adjoint_flux.find(cell_idx);
    if (it_adj == adjoint_flux.end())
      continue;

    double phi_star = it_adj->second;

    // 单能量公式: φ*(r) × χ_p × ν × σ_f × φ(r) × V
    // 使用固定核数据
    sum += phi_star * CHI_PROMPT * NU_TOTAL * SIGMA_F * phi * volume;
  }

  return sum;
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
  std::string nuclear_data_source;
  if (mode_ == BetaEffMode::FIXED_U235) {
    nuclear_data_source = "U-235 thermal (Phase 1 - Fixed)";
  } else {
    nuclear_data_source = "Phase 3 - Dynamic extraction from OpenMC library";
  }
  write_attribute(metadata_group, "nuclear_data_source", nuclear_data_source);
  write_attribute(metadata_group, "energy_groups", n_energy_groups_);
  if (!energy_edges_common_.empty()) {
    write_dataset(metadata_group, "energy_edges", energy_edges_common_);
  }
  write_attribute(metadata_group, "beta_eff_mode",
    mode_ == BetaEffMode::FIXED_U235 ? "FIXED_U235" : "MATERIAL_DEPENDENT");

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

  // 理论值对比 (仅 Phase 1)
  if (mode_ == BetaEffMode::FIXED_U235) {
    double sum_nu_delayed =
      std::accumulate(NU_DELAYED.begin(), NU_DELAYED.end(), 0.0);
    double beta_theoretical = sum_nu_delayed / NU_TOTAL;
    write_attribute(diag_group, "beta_theoretical", beta_theoretical);

    double relative_diff =
      std::abs(beta_total_ - beta_theoretical) / beta_theoretical * 100.0;
    write_attribute(diag_group, "relative_difference_percent", relative_diff);

    // 核数据参数
    write_attribute(diag_group, "nu_total", NU_TOTAL);
    write_attribute(diag_group, "nu_prompt", NU_PROMPT);
    std::vector<double> nu_delayed_vec(NU_DELAYED.begin(), NU_DELAYED.end());
    write_dataset(diag_group, "nu_delayed", nu_delayed_vec);
  } else {
    // Phase 2: 输出材料信息
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
      write_attribute(diag_group, "n_unique_materials",
        static_cast<int>(material_ids.size()));
      write_attribute(diag_group, "total_fissionable_cells",
        static_cast<int>(cell_to_material_.size()));
    }
  }

  H5Gclose(diag_group);
  file_close(file_id);

  std::cout << "  Results written successfully" << std::endl;

  if (mode_ == BetaEffMode::FIXED_U235) {
    double sum_nu_delayed =
      std::accumulate(NU_DELAYED.begin(), NU_DELAYED.end(), 0.0);
    double beta_theoretical = sum_nu_delayed / NU_TOTAL;
    double relative_diff =
      std::abs(beta_total_ - beta_theoretical) / beta_theoretical * 100.0;

    std::cout << "  Theoretical β = " << std::fixed << std::setprecision(5)
              << beta_theoretical << std::endl;
    std::cout << "  Computed β_eff = " << beta_total_ << std::endl;
    std::cout << "  Relative difference = " << std::setprecision(2)
              << relative_diff << "%" << std::endl;
  } else {
    std::cout << "  β_eff = " << std::fixed << std::setprecision(5)
              << beta_total_ << std::endl;
    std::cout << "  Based on " << unique_materials_.size()
              << " fissionable material(s)" << std::endl;
  }
}
//==============================================================================
// Phase 2: Material-Dependent Methods
//==============================================================================

//------------------------------------------------------------------------------
// 材料感知的缓发分子：
//   Σ_cell ∫ φ*(r,E) χ_d,i(E,r) ν_d,i(E,r) Σ_f(r,E) φ(r,E) dE × ΔV。
// 代码重用稀疏网格 map（φ、φ*），若存在单元多群通量则逐群积分，使网格结果
// 能与 Phase 3 提取的材料核数据混合使用；若无多群信息，则退化为基于均匀核
// 数据的单群近似。
//------------------------------------------------------------------------------
double BetaEffective::compute_delayed_numerator_material(int group,
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  const bool use_multi_group =
    n_energy_groups_ > 1 && flux_has_group_data_ && adjoint_has_group_data_ &&
    !flux_group_map_.empty() && !adjoint_group_map_.empty();

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

    auto group_value = [&](const std::vector<double>& values, double fallback,
                         int g) {
      if (values.size() == static_cast<size_t>(n_energy_groups_)) {
        return values[g];
      }
      return fallback;
    };

    auto delayed_value = [&](const std::vector<double>& values,
                           const std::array<double, 8>& fallback, int g) {
      size_t expected =
        static_cast<size_t>(n_energy_groups_) * N_DELAYED_GROUPS;
      if (values.size() == expected) {
        return values[delayed_offset(g, group)];
      }
      return fallback[group];
    };

    for (int g = 0; g < n_energy_groups_; ++g) {
      double phi_g = flux_groups[g];
      double phi_star_g = adjoint_groups[g];
      if (phi_g == 0.0 || phi_star_g == 0.0)
        continue;

      double sigma_f_g =
        group_value(nuc_data.sigma_f_groups, nuc_data.sigma_f, g);
      double nu_delayed_g =
        delayed_value(nuc_data.nu_delayed_groups, nuc_data.nu_delayed, g);
      double chi_delayed_g =
        delayed_value(nuc_data.chi_delayed_groups, nuc_data.chi_delayed, g);

      double term =
        phi_star_g * chi_delayed_g * nu_delayed_g * sigma_f_g * phi_g * volume;
      terms.push_back(term);
    }
  }

  return kahan_sum(terms);
}

//------------------------------------------------------------------------------

// 分母与上式互补：
//   Σ_cell ∫ φ*(r,E) χ_p(E,r) ν_total(E,r) Σ_f(r,E) φ(r,E) dE × ΔV。
// 保持与分子完全对称，从而无论用户是否提供多群网格，β_i = numerator_i /
// denominator 都具备一致的含义。
double BetaEffective::compute_denominator_material(
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  const bool use_multi_group =
    n_energy_groups_ > 1 && flux_has_group_data_ && adjoint_has_group_data_ &&
    !flux_group_map_.empty() && !adjoint_group_map_.empty();

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
        phi_star * nuc_data.nu_total * nuc_data.sigma_f * phi * volume;
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

    auto group_value = [&](const std::vector<double>& values, double fallback,
                         int g) {
      if (values.size() == static_cast<size_t>(n_energy_groups_)) {
        return values[g];
      }
      return fallback;
    };

    for (int g = 0; g < n_energy_groups_; ++g) {
      double phi_g = flux_groups[g];
      double phi_star_g = adjoint_groups[g];
      if (phi_g == 0.0 || phi_star_g == 0.0)
        continue;

      double sigma_f_g =
        group_value(nuc_data.sigma_f_groups, nuc_data.sigma_f, g);
      double nu_total_g =
        group_value(nuc_data.nu_total_groups, nuc_data.nu_total, g);
      double chi_prompt_g =
        group_value(nuc_data.chi_prompt_groups, nuc_data.chi_prompt, g);

      double term =
        phi_star_g * chi_prompt_g * nu_total_g * sigma_f_g * phi_g * volume;
      terms.push_back(term);
    }
  }

  return kahan_sum(terms);
}

//------------------------------------------------------------------------------
// Phase 3: 从 OpenMC 核数据库动态提取核参数
//------------------------------------------------------------------------------
// 通过遍历 OpenMC 的材料与核素表构建位置相关的核数据。对于每个裂变核素，
// 计算宏观 Σ_f、ν_total/ν_prompt/ν_delayed 以及（可选）与公共能群一致的 χ
// 分布。MaterialNuclearData 让后续网格循环无需再次查询物理量即可直接使用。
//------------------------------------------------------------------------------
MaterialNuclearData BetaEffective::extract_material_nuclear_data(
  int material_id) const
{
  MaterialNuclearData data;
  data.material_id = material_id;

  // 安全检查: 确保材料库已初始化
  if (model::materials.empty()) {
    warning("Material library is empty in extract_material_nuclear_data");
    return data;
  }

  // 通过 ID 查找材料 (不是数组索引!)
  Material* mat = nullptr;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->id() == material_id) {
      mat = mat_ptr.get();
      break;
    }
  }

  if (!mat) {
    warning("Material not found for id: " + std::to_string(material_id));
    return data;
  }

  data.material_name = mat->name();
  data.is_fissionable = mat->fissionable();

  if (!data.is_fissionable) {
    return data; // 非裂变材料，返回零值
  }

  const int n_groups = std::max(1, n_energy_groups_);
  const size_t delayed_group_size =
    static_cast<size_t>(n_groups) * N_DELAYED_GROUPS;

  std::vector<double> chi_prompt_acc(n_groups, 0.0);
  std::vector<double> chi_delayed_acc(delayed_group_size, 0.0);
  bool prompt_sampled = false;
  std::array<bool, N_DELAYED_GROUPS> delayed_sampled {};
  delayed_sampled.fill(false);

  data.chi_prompt = 1.0;
  data.chi_delayed.fill(1.0);

  // Phase 3: 使用 OpenMC 核数据库 API 提取能量相关核参数
  // 参考能量与温度由构造函数提供
  const double E_ref = ref_energy_;      // eV
  const double T_ref = ref_temperature_; // K (目前仅用于日志)

  auto sample_spectrum = [&](const ReactionProduct& product,
                           std::vector<double>& spectrum, uint64_t seed_base) {
    if (product.distribution_.empty())
      return false;

    uint64_t seed = seed_base;
    const int n_samples = 5000;
    for (int s = 0; s < n_samples; ++s) {
      double E_out = 0.0;
      double mu = 0.0;
      product.distribution_[0]->sample(E_ref, E_out, mu, &seed);
      int g_idx = group_index_from_energy(E_out);
      if (g_idx >= 0 && g_idx < n_groups) {
        spectrum[g_idx] += 1.0;
      }
    }

    double total = std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
    if (total <= 0.0) {
      double uniform = 1.0 / static_cast<double>(spectrum.size());
      std::fill(spectrum.begin(), spectrum.end(), uniform);
      return false;
    }

    for (double& value : spectrum) {
      value /= total;
    }
    return true;
  };

  // 对材料中的每个裂变核素累加核数据
  // 使用原子密度和裂变截面加权平均: ν̄ = Σ(N_i σ_f,i ν_i) / Σ(N_i σ_f,i)
  double total_fission_density = 0.0; // Σ(N_i × σ_f,i)

  auto nuclides = mat->nuclides();
  auto densities = mat->densities();

  if (nuclides.size() == 0) {
    warning("Material " + std::to_string(material_id) + " has no nuclides");
    return data;
  }

  for (size_t i = 0; i < nuclides.size(); ++i) {
    int nuc_idx = nuclides[i];
    double atom_density = densities[i]; // atom/b-cm

    // 安全检查: 验证核素索引有效性
    if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
      warning("Invalid nuclide index: " + std::to_string(nuc_idx));
      continue;
    }

    // 访问核素数据
    const auto& nuc = data::nuclides[nuc_idx];

    if (!nuc) {
      warning("Nuclide pointer is null for index: " + std::to_string(nuc_idx));
      continue;
    }

    if (!nuc->fissionable_)
      continue;

    // Phase 3: 从 OpenMC 核数据库提取参数
    // 1. 获取裂变截面 σ_f(E, T)
    // 使用 collapse_rate() 方法: 在单能量点计算截面
    // collapse_rate(MT, temperature, energy, flux) -> 返回反应率
    // MT=18 代表裂变反应

    double sigma_f_nuc = 0.0; // 微观裂变截面 (barns)

    try {
      // 方法: 使用 collapse_rate 计算单能群截面
      // 注意: energy.size() 必须等于 flux.size() + 1 (能量边界)
      // 对于单能群，需要 2 个能量边界和 1 个通量权重
      // 使用窄能群近似单能量点: [E_ref - δE, E_ref + δE]
      const double dE = E_ref * 0.001; // 0.1% 能量宽度
      std::vector<double> energy_grid = {E_ref - dE, E_ref + dE}; // 2 个边界
      std::vector<double> flux_weight = {1.0}; // 1 个能群的通量权重

      // collapse_rate 参数: (MT, temperature [K], energy [eV], flux)
      const int MT_FISSION = 18;

      // 直接使用参考温度 T_ref (K)
      sigma_f_nuc =
        nuc->collapse_rate(MT_FISSION, T_ref, energy_grid, flux_weight);

    } catch (...) {
      // 如果 collapse_rate 失败，尝试备用方法
      warning("Failed to extract cross section for " + nuc->name_ +
              ", using fallback method");

      // 备用方法: 直接查表 (需要能量网格索引)
      // 这里使用热中子典型值作为后备
      std::string nuc_name = nuc->name_;
      if (nuc_name.find("U235") != std::string::npos) {
        sigma_f_nuc = 584.4; // barns at 0.0253 eV
      } else if (nuc_name.find("Pu239") != std::string::npos) {
        sigma_f_nuc = 747.4;
      } else if (nuc_name.find("Pu241") != std::string::npos) {
        sigma_f_nuc = 1009.0;
      } else {
        sigma_f_nuc = 584.4; // 默认使用 U-235
      }
    }

    // 2. 获取中子产额 ν(E)
    // 使用 Nuclide::nu(E, mode, group) 方法
    // EmissionMode: prompt (瞬发), delayed (缓发), total (总和)

    double nu_total_nuc = 0.0;
    double nu_prompt_nuc = 0.0;
    std::array<double, 8> nu_delayed_nuc = {}; // Phase 3: 8组

    try {
      // 获取总中子产额
      nu_total_nuc = nuc->nu(E_ref, ReactionProduct::EmissionMode::total, 0);

      // 获取瞬发中子产额
      nu_prompt_nuc = nuc->nu(E_ref, ReactionProduct::EmissionMode::prompt, 0);

      // 获取缓发中子产额
      // Phase 3 改进: 直接访问 products_ 数组，检查粒子类型和发射模式
      // products_[0] = 瞬发中子
      // products_[1..N] = 可能是缓发中子、光子等其他产物

      // 先获取总缓发产额用于验证
      double nu_delayed_total =
        nuc->nu(E_ref, ReactionProduct::EmissionMode::delayed, 0);

      if (nu_delayed_total < 1e-10) {
        // 没有缓发中子数据
        nu_delayed_nuc.fill(0.0);
      } else {
        // 遍历裂变反应的产物，只提取缓发中子
        if (!nuc->fission_rx_.empty()) {
          const auto& fission = nuc->fission_rx_[0];
          int delayed_group_idx = 0; // 缓发群索引 (0-7)

          for (int i = 1;
            i < fission->products_.size() && delayed_group_idx < 8; ++i) {
            const auto& product = fission->products_[i];

            // 检查是否是中子且是缓发发射
            if (product.particle_ == ParticleType::neutron &&
                product.emission_mode_ ==
                  ReactionProduct::EmissionMode::delayed) {

              // 提取该组的中子产额
              double nu_d = (*product.yield_)(E_ref);
              nu_delayed_nuc[delayed_group_idx] = nu_d;
              delayed_group_idx++;
            }
          }

          // 剩余组设为0
          for (int g = delayed_group_idx; g < 8; ++g) {
            nu_delayed_nuc[g] = 0.0;
          }
        } else {
          // 没有裂变反应数据
          nu_delayed_nuc.fill(0.0);
        }
      }

    } catch (...) {
      // 如果 nu() 方法失败，使用备用硬编码数据
      warning(
        "Failed to extract nu for " + nuc->name_ + ", using fallback data");

      std::string nuc_name = nuc->name_;
      if (nuc_name.find("U235") != std::string::npos) {
        nu_total_nuc = 2.43;
        nu_prompt_nuc = 2.42;
        // U-235 典型6组数据，第7-8组为0
        nu_delayed_nuc = {
          0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273, 0.0, 0.0};
      } else if (nuc_name.find("Pu239") != std::string::npos) {
        nu_total_nuc = 2.88;
        nu_prompt_nuc = 2.86;
        nu_delayed_nuc = {
          0.000065, 0.000539, 0.000425, 0.000799, 0.000277, 0.000091, 0.0, 0.0};
      } else if (nuc_name.find("Pu241") != std::string::npos) {
        nu_total_nuc = 2.93;
        nu_prompt_nuc = 2.91;
        nu_delayed_nuc = {
          0.000066, 0.000496, 0.000471, 0.000869, 0.000259, 0.000079, 0.0, 0.0};
      } else {
        // 默认使用 U-235
        nu_total_nuc = 2.43;
        nu_prompt_nuc = 2.42;
        nu_delayed_nuc = {
          0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273, 0.0, 0.0};
      }
    }

    // Phase 4: 提取中子能谱权重 —— 单能群版本暂时停用
    // 原因: 在单群公式里直接乘 χ(E) 概率密度会导致维度不自洽
    // 多群实现时会使用群常数 χ_g，而非连续能谱的 χ(E_ref)

    // ★ 单能群版本：不提取 χ，避免重复能谱折合
    /*
    double chi_prompt_nuc = 1.0;
    std::array<double, 8> chi_delayed_nuc = {};
    chi_delayed_nuc.fill(1.0);

    try {
      if (!nuc->fission_rx_.empty()) {
        const auto& fission = nuc->fission_rx_[0];

        const int n_samples = 10000;     // 每个能谱的总采样数
        const double E_incident = 1.0e6; // 1 MeV 入射能量

        // 瞬发中子参考能量: 2.0 MeV (裂变中子特征能量)
        const double E_ref_prompt = 2.0e6;           // eV
        const double dE_prompt = E_ref_prompt * 0.5; // ±50% 窗口
        const double E_low_prompt = E_ref_prompt - dE_prompt;
        const double E_high_prompt = E_ref_prompt + dE_prompt;

        // 1. 采样瞬发中子能谱
        for (const auto& product : fission->products_) {
          if (product.particle_ == ParticleType::neutron &&
              product.emission_mode_ == ReactionProduct::EmissionMode::prompt) {

            if (!product.distribution_.empty()) {
              uint64_t seed = 987654321ULL;
              int count_in_window = 0;

              for (int s = 0; s < n_samples; ++s) {
                double E_out = 0.0;
                double mu = 0.0;
                product.distribution_[0]->sample(E_incident, E_out, mu, &seed);

                if (E_out >= E_low_prompt && E_out <= E_high_prompt) {
                  count_in_window++;
                }
              }

              // 概率密度 = 命中数 / (总采样数 × 窗口宽度)
              if (count_in_window > 0) {
                chi_prompt_nuc = static_cast<double>(count_in_window) /
                                 (n_samples * 2.0 * dE_prompt);
              } else {
                // 使用 Maxwell 分布解析估计: χ(E) ∝ √E exp(-E/T)
                // T_prompt ≈ 1.29 MeV
                const double T_prompt = 1.29e6;
                chi_prompt_nuc = 0.484 * std::sqrt(E_ref_prompt) *
                                 std::exp(-E_ref_prompt / T_prompt);
              }
            }
            break;
          }
        }

        // 缓发中子参考能量: 0.5 MeV (缓发中子特征能量)
        const double E_ref_delayed = 0.5e6;            // eV
        const double dE_delayed = E_ref_delayed * 0.5; // ±50% 窗口
        const double E_low_delayed = E_ref_delayed - dE_delayed;
        const double E_high_delayed = E_ref_delayed + dE_delayed;

        // 2. 采样缓发中子能谱
        std::vector<const ReactionProduct*> delayed_products;
        delayed_products.reserve(8);

        for (const auto& product : fission->products_) {
          if (product.particle_ == ParticleType::neutron &&
              product.emission_mode_ ==
                ReactionProduct::EmissionMode::delayed) {
            delayed_products.push_back(&product);
          }
        }

        if (!delayed_products.empty() && delayed_products.size() <= 8) {
          uint64_t seed = 123456789ULL + i;

          for (size_t g = 0; g < delayed_products.size(); ++g) {
            const auto* product = delayed_products[g];

            if (product->distribution_.empty())
              continue;

            int count_in_window = 0;

            for (int s = 0; s < n_samples; ++s) {
              double E_out = 0.0;
              double mu = 0.0;
              product->distribution_[0]->sample(E_incident, E_out, mu, &seed);

              if (E_out >= E_low_delayed && E_out <= E_high_delayed) {
                count_in_window++;
              }
            }

            // 概率密度 = 命中数 / (总采样数 × 窗口宽度)
            if (count_in_window > 0) {
              chi_delayed_nuc[g] = static_cast<double>(count_in_window) /
                                   (n_samples * 2.0 * dE_delayed);
            } else {
              // 使用 Maxwell 分布解析估计: χ(E) ∝ √E exp(-E/T)
              // T_delayed ≈ 0.4 MeV (缓发中子温度更低)
              const double T_delayed = 0.4e6;
              chi_delayed_nuc[g] = 0.484 * std::sqrt(E_ref_delayed) *
                                   std::exp(-E_ref_delayed / T_delayed);
            }
          }
        }
      }
    } catch (...) {
      warning("Failed to extract chi spectrum for " + nuc->name_ +
              ", using default values");
      chi_prompt_nuc = 1.0;
      chi_delayed_nuc.fill(1.0);
    }
    */

    // 3. 累加加权核数据
    // 宏观裂变截面贡献: N_i × σ_f,i
    // atom_density 单位: atom/b-cm (OpenMC 标准单位)
    // sigma_f_nuc 单位: barn
    // 结果单位: (atom/b-cm) × barn = atom/cm = cm⁻¹ (宏观截面)
    double macro_sigma_f = atom_density * sigma_f_nuc; // 不需要额外转换！

    if (n_groups > 1 && !energy_edges_common_.empty() &&
        !nuc->fission_rx_.empty()) {
      const auto& fission = nuc->fission_rx_[0];

      // 采样瞬发能谱
      for (const auto& product : fission->products_) {
        if (product.particle_ == ParticleType::neutron &&
            product.emission_mode_ == ReactionProduct::EmissionMode::prompt) {
          std::vector<double> prompt_spectrum(n_groups, 0.0);
          uint64_t seed = 0x9e3779b97f4a7c15ULL ^
                          static_cast<uint64_t>(material_id) ^
                          (static_cast<uint64_t>(nuc_idx) << 16);
          sample_spectrum(product, prompt_spectrum, seed);
          for (int g = 0; g < n_groups; ++g) {
            chi_prompt_acc[g] += macro_sigma_f * prompt_spectrum[g];
          }
          prompt_sampled = true;
          break;
        }
      }

      // 采样缓发能谱
      int delayed_group_idx = 0;
      for (const auto& product : fission->products_) {
        if (product.particle_ == ParticleType::neutron &&
            product.emission_mode_ == ReactionProduct::EmissionMode::delayed) {
          std::vector<double> delayed_spectrum(n_groups, 0.0);
          uint64_t seed = 0x6a09e667f3bcc909ULL ^
                          static_cast<uint64_t>(material_id) ^
                          (static_cast<uint64_t>(nuc_idx) << 24) ^
                          static_cast<uint64_t>(delayed_group_idx);
          sample_spectrum(product, delayed_spectrum, seed);
          for (int g = 0; g < n_groups; ++g) {
            chi_delayed_acc[delayed_offset(g, delayed_group_idx)] +=
              macro_sigma_f * delayed_spectrum[g];
          }
          delayed_sampled[delayed_group_idx] = true;
          delayed_group_idx++;
          if (delayed_group_idx >= N_DELAYED_GROUPS)
            break;
        }
      }
    }

    total_fission_density += macro_sigma_f;

    data.sigma_f += macro_sigma_f;
    data.nu_total += macro_sigma_f * nu_total_nuc;
    data.nu_prompt += macro_sigma_f * nu_prompt_nuc;

    for (int g = 0; g < 8; ++g) { // Phase 3: 8组
      data.nu_delayed[g] += macro_sigma_f * nu_delayed_nuc[g];
      // ★ 单能群版本：不累积 chi_delayed
      // data.chi_delayed[g] += macro_sigma_f * chi_delayed_nuc[g];
    }

    // ★ 单能群版本：不累积 chi_prompt
    // data.chi_prompt += macro_sigma_f * chi_prompt_nuc;
  }

  // 4. 归一化 (除以总裂变密度)
  if (total_fission_density > 0.0) {
    data.nu_total /= total_fission_density;
    data.nu_prompt /= total_fission_density;

    // 归一化中子产额
    for (int g = 0; g < 8; ++g) { // Phase 3: 8组
      data.nu_delayed[g] /= total_fission_density;
    }

    if (prompt_sampled) {
      for (double& value : chi_prompt_acc) {
        value /= total_fission_density;
      }
    }
    if (!chi_delayed_acc.empty()) {
      for (double& value : chi_delayed_acc) {
        value /= total_fission_density;
      }
    }

    // ★ 单能群版本：不归一化 chi_prompt / chi_delayed
    // data.chi_prompt /= total_fission_density;
    // for (int g = 0; g < 8; ++g) {
    //   data.chi_delayed[g] /= total_fission_density;
    // }

    // 能谱按裂变密度加权，保持 χ_g/χ_{d,ig} 的概率意义

    // sigma_f 已经是宏观截面，不需要归一化
  }

  if (!prompt_sampled || total_fission_density <= 0.0) {
    double uniform = 1.0 / static_cast<double>(n_groups);
    std::fill(chi_prompt_acc.begin(), chi_prompt_acc.end(), uniform);
  }

  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    if (!delayed_sampled[d] || total_fission_density <= 0.0) {
      double uniform = 1.0 / static_cast<double>(n_groups);
      for (int g = 0; g < n_groups; ++g) {
        chi_delayed_acc[delayed_offset(g, d)] = uniform;
      }
    }
  }

  data.chi_prompt_groups = chi_prompt_acc;
  data.chi_delayed_groups = chi_delayed_acc;

  data.sigma_f_groups.assign(n_groups, data.sigma_f);
  data.nu_total_groups.assign(n_groups, data.nu_total);
  data.nu_prompt_groups.assign(n_groups, data.nu_prompt);

  data.nu_delayed_groups.assign(delayed_group_size, 0.0);
  for (int g = 0; g < n_groups; ++g) {
    for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
      data.nu_delayed_groups[delayed_offset(g, d)] = data.nu_delayed[d];
    }
  }

  return data;
}

//------------------------------------------------------------------------------

void BetaEffective::build_cell_material_map(
  const std::unordered_map<int, double>& flux)
{
  std::cout << "  Building cell-material mapping using geometry queries..."
            << std::endl;

  cell_to_material_.clear();
  cell_nuclear_data_.clear();
  n_heterogeneous_cells_ = 0;

  // 检查材料库是否可用
  if (model::materials.empty()) {
    warning("Material library is empty! Cannot perform geometry query.");
    return;
  }

  std::cout << "  Total materials in library: " << model::materials.size()
            << std::endl;
  std::cout << "  Querying geometry for " << flux.size() << " mesh cells..."
            << std::endl;

  // 线程局部缓存和统计
  struct ThreadLocalData {
    std::unordered_map<int, MaterialNuclearData> material_data_cache;
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

#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    // 初始化线程局部存储
#pragma omp single
    {
      thread_data.resize(num_threads);
    }

    auto& local_data = thread_data[thread_id];
    GeometryState geom; // 线程私有几何状态

    // 将 flux map 转换为 vector 以便并行化 (读取顺序与稀疏数据一致)
    std::vector<std::pair<int, double>> flux_vec(flux.begin(), flux.end());

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

      // 对每个网格单元执行 1/8/27 点采样来估计材料体积分数（Phase 2.3 核心）。
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
        // 同质单元: 直接使用单一材料的数据
        int mat_id = material_counts.begin()->first;

        // 缓存材料核数据
        if (local_data.material_data_cache.find(mat_id) ==
            local_data.material_data_cache.end()) {
          local_data.material_data_cache[mat_id] =
            extract_material_nuclear_data(mat_id);
        }

        cell_data = local_data.material_data_cache[mat_id];
        local_data.cell_to_material[cell_idx] = mat_id;
      } else {
        // 异质单元: 计算加权核数据(仅基于裂变材料)
        cell_data =
          compute_weighted_nuclear_data(material_counts, fissionable_hits);
        local_data.cell_to_material[cell_idx] = -1; // 标记为多材料

        // 缓存涉及的材料数据
        for (const auto& [mat_id, count] : material_counts) {
          if (local_data.material_data_cache.find(mat_id) ==
              local_data.material_data_cache.end()) {
            local_data.material_data_cache[mat_id] =
              extract_material_nuclear_data(mat_id);
          }
        }
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
  std::unordered_map<int, MaterialNuclearData> material_data_cache;
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

    // 合并材料缓存
    for (const auto& [mat_id, data] : local.material_data_cache) {
      material_data_cache[mat_id] = data;
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

  // 构建唯一材料列表
  unique_materials_.clear();
  for (const auto& [mat_id, data] : material_data_cache) {
    unique_materials_.push_back(data);
  }

  // 输出统计信息
  std::cout << "  " << std::string(60, '-') << std::endl;
  std::cout << "  Geometry query results:" << std::endl;
  std::cout << "    Fissionable cells: " << fissionable_count << std::endl;
  if (n_sample_points_ > 1) {
    std::cout << "    Heterogeneous cells (multi-fissionable): "
              << n_heterogeneous_cells_ << std::endl;
  }
  std::cout << "    Non-fissionable cells: " << non_fissionable_count
            << std::endl;
  std::cout << "    Void cells: " << void_count << std::endl;
  if (geometry_failed_count > 0) {
    std::cout << "    Geometry query failed: " << geometry_failed_count
              << std::endl;
  }
  std::cout << "    Unique fissionable materials: " << unique_materials_.size()
            << std::endl;

  // 输出每个裂变材料的分布(仅计算裂变区域采样)
  if (!material_hit_counts.empty()) {
    std::cout << "\n  Fissionable material distribution (sample hits):"
              << std::endl;
    for (const auto& [mat_id, count] : material_hit_counts) {
      if (material_data_cache.find(mat_id) != material_data_cache.end()) {
        const auto& data = material_data_cache[mat_id];
        std::cout << "    Material " << mat_id << " (" << data.material_name
                  << "): " << count << " hits" << std::endl;
        std::cout << "      ν_total = " << std::fixed << std::setprecision(4)
                  << data.nu_total << ", Σ_f = " << std::scientific
                  << std::setprecision(3) << data.sigma_f << " cm⁻¹"
                  << std::endl;

        // ★ 单能群版本：不再输出 χ 信息
        // 多群版本实现后会输出每群的 χ_p,g 和 χ_d,ig
      }
    }
  }
  std::cout << "  " << std::string(60, '-') << std::endl;
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
  const std::unordered_map<int, int>& material_counts, int total_samples) const
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

      // 提取材料核数据
      MaterialNuclearData mat_data = extract_material_nuclear_data(mat_id);

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

      // 提取材料核数据
      MaterialNuclearData mat_data = extract_material_nuclear_data(mat_id);

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
  // 参考: Kahan, W. (1965). "Further remarks on reducing truncation errors"
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