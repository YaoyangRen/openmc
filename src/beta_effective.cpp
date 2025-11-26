#include "openmc/beta_effective.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"

namespace openmc {

void BetaEffective::compute_from_files(const std::string& flux_file,
  const std::string& adjoint_flux_file, const std::string& output_file)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "EFFECTIVE DELAYED NEUTRON FRACTION COMPUTATION" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

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

  std::cout << "\n[3/4] Computing β_eff using U-235 nuclear data" << std::endl;
  std::cout << "  ν_total = " << NU_TOTAL << std::endl;
  std::cout << "  ν_prompt = " << NU_PROMPT << std::endl;
  std::cout << "  Cell volume = " << volume << " cm³" << std::endl;

  // 4. 计算分母(只需计算一次)
  denominator_ = compute_denominator(flux, adjoint_flux, volume);

  if (denominator_ <= 0.0) {
    fatal_error("Denominator is zero or negative! Cannot compute β_eff.");
  }

  std::cout << "  Denominator = " << std::scientific << std::setprecision(6)
            << denominator_ << std::endl;

  // 5. 对每个缓发群计算分子和 β_i,eff
  std::cout << "\n  Delayed group contributions:" << std::endl;
  std::cout << "  Group    ν_d,i      Numerator        β_i,eff" << std::endl;
  std::cout << "  " << std::string(55, '-') << std::endl;

  for (int i = 0; i < 6; ++i) {
    numerators_[i] = compute_delayed_numerator(i, flux, adjoint_flux, volume);
    beta_i_[i] = numerators_[i] / denominator_;

    std::cout << "    " << (i + 1) << "    " << std::scientific
              << std::setprecision(6) << NU_DELAYED[i] << "   "
              << numerators_[i] << "   " << beta_i_[i] << std::endl;
  }

  // 6. 计算总 β_eff
  beta_total_ = std::accumulate(beta_i_.begin(), beta_i_.end(), 0.0);

  std::cout << "  " << std::string(55, '-') << std::endl;
  std::cout << "  Total β_eff = " << std::fixed << std::setprecision(5)
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
  double& pitch) const
{
  hid_t file_id = file_open(filename, 'r');

  // 读取网格参数
  std::array<int, 3> grid_shape;
  read_dataset(file_id, "grid_shape", grid_shape);
  shape = grid_shape;

  std::array<double, 3> grid_pitch_array;
  read_dataset(file_id, "grid_pitch", grid_pitch_array);
  pitch = grid_pitch_array[0]; // 假设各向同性

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
}

//------------------------------------------------------------------------------

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
  // beta_i [6] - 各组 β_i,eff
  std::vector<double> beta_i_vec(beta_i_.begin(), beta_i_.end());
  write_dataset(file_id, "beta_i", beta_i_vec);

  // beta_total - 总 β_eff
  std::vector<double> beta_total_vec = {beta_total_};
  write_dataset(file_id, "beta_total", beta_total_vec);

  // uncertainty [6] - Phase 1暂时为零,预留给Phase 4统计不确定度
  std::vector<double> uncertainty(6, 0.0);
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
  write_attribute(
    metadata_group, "nuclear_data_source", "U-235 thermal (Phase 1)");
  write_attribute(metadata_group, "energy_groups", 1); // 单能群
  write_dataset(metadata_group, "grid_shape", grid_shape_);
  write_attribute(metadata_group, "grid_pitch", grid_pitch_);

  H5Gclose(metadata_group);

  // ========== diagnostics/ 组 ==========
  hid_t diag_group =
    H5Gcreate(file_id, "diagnostics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // numerator [6] - 各组分子项
  std::vector<double> numerators_vec(numerators_.begin(), numerators_.end());
  write_dataset(diag_group, "numerator", numerators_vec);

  // denominator - 分母项(单个值)
  std::vector<double> denom_vec = {denominator_};
  write_dataset(diag_group, "denominator", denom_vec);

  // normalization_factor - 归一化因子(体积元)
  write_attribute(diag_group, "normalization_factor",
    grid_pitch_ * grid_pitch_ * grid_pitch_);

  // 理论值对比
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

  H5Gclose(diag_group);
  file_close(file_id);

  std::cout << "  Results written successfully" << std::endl;
  std::cout << "  Theoretical β = " << std::fixed << std::setprecision(5)
            << beta_theoretical << std::endl;
  std::cout << "  Computed β_eff = " << beta_total_ << std::endl;
  std::cout << "  Relative difference = " << std::setprecision(2)
            << relative_diff << "%" << std::endl;
}

} // namespace openmc
