#include "openmc/beta_effective.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/reaction_product.h" // Phase 3: EmissionMode 枚举
#include "openmc/settings.h"

namespace openmc {

void BetaEffective::compute_from_files(const std::string& flux_file,
  const std::string& adjoint_flux_file, const std::string& output_file)
{
  std::cout << std::string(70, '=') << std::endl;
  std::cout << "EFFECTIVE DELAYED NEUTRON FRACTION COMPUTATION" << std::endl;

  // 1. 读取正向通量数据
  std::cout << "\n[1/4] Reading forward flux from: " << flux_file << std::endl;

  if (!file_exists(flux_file)) {
    fatal_error("Flux file not found: " + flux_file +
                "\nPlease ensure flux_mesh has been computed.");
  }

  std::unordered_map<int, double> flux;
  std::array<int, 3> flux_shape;
  double flux_pitch;
  read_flux_data(flux_file, flux, flux_shape, flux_pitch);

  std::cout << "  Grid: " << flux_shape[0] << " x " << flux_shape[1] << " x "
            << flux_shape[2] << std::endl;
  std::cout << "  Pitch: " << flux_pitch << " cm" << std::endl;
  std::cout << "  Non-zero cells: " << flux.size() << std::endl;

  // 2. 读取共轭通量数据
  std::cout << "\n[2/4] Reading adjoint flux from: " << adjoint_flux_file
            << std::endl;

  if (!file_exists(adjoint_flux_file)) {
    fatal_error("Adjoint flux file not found: " + adjoint_flux_file +
                "\nPlease ensure adjoint flux has been computed.");
  }

  std::unordered_map<int, double> adjoint_flux;
  std::array<int, 3> adjoint_shape;
  double adjoint_pitch;
  read_adjoint_flux_data(
    adjoint_flux_file, adjoint_flux, adjoint_shape, adjoint_pitch);

  std::cout << "  Grid: " << adjoint_shape[0] << " x " << adjoint_shape[1]
            << " x " << adjoint_shape[2] << std::endl;
  std::cout << "  Pitch: " << adjoint_pitch << " cm" << std::endl;
  std::cout << "  Non-zero cells: " << adjoint_flux.size() << std::endl;

  // 3. 验证网格一致性
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

  // 计算网格单元体积
  double volume = flux_pitch * flux_pitch * flux_pitch;

  // 根据模式选择计算方法
  if (mode_ == BetaEffMode::FIXED_U235) {
    // ===== Phase 1: 固定 U-235 核数据 =====
    std::cout << "\n[3/4] Computing β_eff using FIXED U-235 nuclear data"
              << std::endl;
    std::cout << "  Mode: Phase 1 (FIXED_U235)" << std::endl;
    std::cout << "  ν_total = " << NU_TOTAL << std::endl;
    std::cout << "  ν_prompt = " << NU_PROMPT << std::endl;
    std::cout << "  Cell volume = " << volume << " cm³" << std::endl;

    // 计算分母
    denominator_ = compute_denominator(flux, adjoint_flux, volume);

    if (denominator_ <= 0.0) {
      fatal_error("Denominator is zero or negative! Cannot compute β_eff.");
    }

    std::cout << "  Denominator = " << std::scientific << std::setprecision(6)
              << denominator_ << std::endl;

    // 对每个缓发群计算分子和 β_i,eff (Phase 3: 8组)
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

    // 构建单元-材料映射并提取核数据
    build_cell_material_map(flux);

    // 计算分母
    denominator_ = compute_denominator_material(flux, adjoint_flux, volume);

    if (denominator_ <= 0.0) {
      fatal_error("Denominator is zero or negative! Cannot compute β_eff.");
    }

    std::cout << "  Denominator = " << std::scientific << std::setprecision(6)
              << denominator_ << std::endl;

    // 对每个缓发群计算分子和 β_i,eff
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

void BetaEffective::read_flux_data(const std::string& filename,
  std::unordered_map<int, double>& flux_map, std::array<int, 3>& shape,
  double& pitch)
{
  hid_t file_id = file_open(filename, 'r');

  // 读取网格参数
  std::array<int, 3> grid_shape;
  read_dataset(file_id, "grid_shape", grid_shape);
  shape = grid_shape;

  std::array<double, 3> grid_pitch_array;
  read_dataset(file_id, "grid_pitch", grid_pitch_array);
  pitch = grid_pitch_array[0]; // 假设各向同性

  // 读取网格左下角坐标 (Phase 2 需要)
  if (mode_ == BetaEffMode::MATERIAL_DEPENDENT) {
    read_dataset(file_id, "grid_lower_left", grid_lower_left_);
  }

  // 读取稀疏通量数据
  std::vector<int> cell_indices;
  std::vector<double> flux_mean;

  read_dataset(file_id, "cell_indices", cell_indices);
  read_dataset(file_id, "flux_mean", flux_mean);

  // 构建稀疏 map
  flux_map.clear();
  for (size_t i = 0; i < cell_indices.size(); ++i) {
    flux_map[cell_indices[i]] = flux_mean[i];
  }

  file_close(file_id);
} //------------------------------------------------------------------------------

void BetaEffective::read_adjoint_flux_data(const std::string& filename,
  std::unordered_map<int, double>& adjoint_flux_map, std::array<int, 3>& shape,
  double& pitch) const
{
  hid_t file_id = file_open(filename, 'r');

  // 读取网格参数
  std::array<int, 3> grid_shape;
  read_dataset(file_id, "shape", grid_shape);
  shape = grid_shape;

  read_attribute(file_id, "pitch", pitch);

  // 读取稀疏共轭通量数据
  std::vector<int> cell_indices;
  std::vector<double> adjoint_flux_values;

  read_dataset(file_id, "cell_indices", cell_indices);
  read_dataset(file_id, "adjoint_flux_values", adjoint_flux_values);

  // 构建稀疏 map
  adjoint_flux_map.clear();
  for (size_t i = 0; i < cell_indices.size(); ++i) {
    adjoint_flux_map[cell_indices[i]] = adjoint_flux_values[i];
  }

  file_close(file_id);
}

//------------------------------------------------------------------------------

double BetaEffective::get_beta_i(int group) const
{
  if (group < 0 || group >= 6) {
    fatal_error("Invalid delayed group index: " + std::to_string(group) +
                " (must be 0-5)");
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
  write_attribute(metadata_group, "energy_groups", 1); // 单能群
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

double BetaEffective::compute_delayed_numerator_material(int group,
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  // 收集所有贡献项用于 Kahan 求和
  std::vector<double> terms;
  terms.reserve(flux.size());

  // 遍历所有有正向通量的单元
  for (const auto& [cell_idx, phi] : flux) {
    // 检查该单元是否也有共轭通量
    auto it_adj = adjoint_flux.find(cell_idx);
    if (it_adj == adjoint_flux.end())
      continue;

    double phi_star = it_adj->second;

    // 获取该单元的核数据
    auto it_data = cell_nuclear_data_.find(cell_idx);
    if (it_data == cell_nuclear_data_.end() || !it_data->second.is_fissionable)
      continue;

    const auto& nuc_data = it_data->second;

    // 使用材料相关核数据: φ*(r) × χ_d,i × ν_d,i(r) × σ_f(r) × φ(r) × V
    double term = phi_star * nuc_data.chi_delayed[group] *
                  nuc_data.nu_delayed[group] * nuc_data.sigma_f * phi * volume;
    terms.push_back(term);
  }

  // 使用 Kahan 求和提高精度
  return kahan_sum(terms);
}

//------------------------------------------------------------------------------

double BetaEffective::compute_denominator_material(
  const std::unordered_map<int, double>& flux,
  const std::unordered_map<int, double>& adjoint_flux, double volume) const
{
  // 收集所有贡献项用于 Kahan 求和
  std::vector<double> terms;
  terms.reserve(flux.size());

  // 遍历所有有正向通量的单元
  for (const auto& [cell_idx, phi] : flux) {
    // 检查该单元是否也有共轭通量
    auto it_adj = adjoint_flux.find(cell_idx);
    if (it_adj == adjoint_flux.end())
      continue;

    double phi_star = it_adj->second;

    // 获取该单元的核数据
    auto it_data = cell_nuclear_data_.find(cell_idx);
    if (it_data == cell_nuclear_data_.end() || !it_data->second.is_fissionable)
      continue;

    const auto& nuc_data = it_data->second;

    // 使用材料相关核数据: φ*(r) × χ_p × ν(r) × σ_f(r) × φ(r) × V
    double term = phi_star * nuc_data.chi_prompt * nuc_data.nu_total *
                  nuc_data.sigma_f * phi * volume;
    terms.push_back(term);
  }

  // 使用 Kahan 求和提高精度
  return kahan_sum(terms);
}

//------------------------------------------------------------------------------
// Phase 3: 从 OpenMC 核数据库动态提取核参数
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

  // Phase 3: 使用 OpenMC 核数据库 API 提取能量相关核参数
  // 参考能量: 0.0253 eV (热中子), 参考温度: 293.6 K
  const double E_ref = ref_energy_;      // eV (构造函数传入)
  const double T_ref = ref_temperature_; // K (构造函数传入)

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

    // 3. 累加加权核数据
    // 宏观裂变截面贡献: N_i × σ_f,i
    double macro_sigma_f = atom_density * sigma_f_nuc * 1e-24; // barn to cm^2

    total_fission_density += macro_sigma_f;

    data.sigma_f += macro_sigma_f;
    data.nu_total += macro_sigma_f * nu_total_nuc;
    data.nu_prompt += macro_sigma_f * nu_prompt_nuc;

    for (int g = 0; g < 8; ++g) { // Phase 3: 8组
      data.nu_delayed[g] += macro_sigma_f * nu_delayed_nuc[g];
    }
  }

  // 4. 归一化 (除以总裂变密度)
  if (total_fission_density > 0.0) {
    data.nu_total /= total_fission_density;
    data.nu_prompt /= total_fission_density;
    for (int g = 0; g < 8; ++g) { // Phase 3: 8组
      data.nu_delayed[g] /= total_fission_density;
    }
    // sigma_f 已经是宏观截面，不需要归一化
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

    // 将 flux map 转换为 vector 以便并行化
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
      int valid_samples = 0;

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
          continue; // 跳过 void
        }

        // 验证材料索引
        if (material_idx < 0 ||
            material_idx >= static_cast<int>(model::materials.size())) {
          continue;
        }

        const auto& material = model::materials[material_idx];
        if (!material || !material->fissionable()) {
          continue; // 跳过非裂变材料
        }

        int material_id = material->id();
        material_counts[material_id]++;
        valid_samples++;
      }

      // 处理采样结果
      if (valid_samples == 0) {
        // 没有有效采样点
        if (material_counts.empty()) {
          local_data.void_count++;
        } else {
          local_data.non_fissionable_count++;
        }
        continue;
      }

      // 检查是否异质
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
        // 异质单元: 计算加权核数据
        cell_data =
          compute_weighted_nuclear_data(material_counts, valid_samples);
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
    std::cout << "    Heterogeneous cells (multi-material): "
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

  // 输出每个材料的分布
  if (!material_hit_counts.empty()) {
    std::cout << "\n  Material distribution (sample hits):" << std::endl;
    for (const auto& [mat_id, count] : material_hit_counts) {
      if (material_data_cache.find(mat_id) != material_data_cache.end()) {
        const auto& data = material_data_cache[mat_id];
        std::cout << "    Material " << mat_id << " (" << data.material_name
                  << "): " << count << " hits" << std::endl;
        std::cout << "      ν_total = " << std::fixed << std::setprecision(4)
                  << data.nu_total << ", Σ_f = " << std::scientific
                  << std::setprecision(3) << data.sigma_f << " cm⁻¹"
                  << std::endl;
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
  MaterialNuclearData weighted_data;
  weighted_data.material_id = -1; // 标记为多材料单元
  weighted_data.material_name = "mixed";
  weighted_data.is_fissionable = true;

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