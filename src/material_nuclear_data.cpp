#include "openmc/material_nuclear_data.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/reaction_product.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>

namespace openmc {

//==============================================================================
// MaterialNuclearDataExtractor 实现
//==============================================================================

MaterialNuclearDataExtractor::MaterialNuclearDataExtractor(int n_energy_groups,
  const std::vector<double>& energy_edges,
  const std::vector<double>& flux_weights, double ref_energy,
  double ref_temperature)
  : n_groups_(n_energy_groups), energy_edges_(energy_edges),
    flux_weights_(flux_weights), ref_energy_(ref_energy),
    ref_temperature_(ref_temperature)
{
  // 校验输入
  if (n_groups_ <= 1) {
    fatal_error("MaterialNuclearDataExtractor requires n_groups > 1");
  }
  if (energy_edges_.size() != static_cast<size_t>(n_groups_ + 1)) {
    fatal_error("energy_edges size must equal n_groups + 1");
  }
  // 若未提供 flux_weights，使用均匀分布
  if (flux_weights_.empty()) {
    flux_weights_.assign(n_groups_, 1.0 / n_groups_);
  }
}

void derive_diagnostic_nuclear_data_from_source_terms(MaterialNuclearData& data,
  int n_energy_groups, const std::vector<double>& collapse_weights)
{
  const int n_groups = n_energy_groups;
  if (n_groups <= 0)
    return;

  const size_t prompt_size =
    static_cast<size_t>(n_groups) * static_cast<size_t>(n_groups);
  const size_t delayed_size = static_cast<size_t>(n_groups) *
                              N_DELAYED_GROUPS *
                              static_cast<size_t>(n_groups);
  const size_t delayed_legacy_size =
    static_cast<size_t>(n_groups) * N_DELAYED_GROUPS;

  if (data.prompt_prod_groups.size() != prompt_size ||
      data.delayed_prod_groups.size() != delayed_size) {
    return;
  }

  if (data.sigma_f_groups.size() != static_cast<size_t>(n_groups)) {
    data.sigma_f_groups.assign(n_groups, 0.0);
  }

  data.nu_total_groups.assign(n_groups, 0.0);
  data.nu_prompt_groups.assign(n_groups, 0.0);
  data.nu_delayed_groups.assign(delayed_legacy_size, 0.0);
  data.chi_prompt_groups.assign(n_groups, 0.0);
  data.chi_delayed_groups.assign(delayed_legacy_size, 0.0);
  data.nu_delayed.fill(0.0);
  data.chi_delayed.fill(0.0);

  std::vector<double> prompt_source_by_group(n_groups, 0.0);
  std::vector<double> total_source_by_group(n_groups, 0.0);
  std::vector<double> delayed_source_by_group(delayed_legacy_size, 0.0);

  double total_prompt_source = 0.0;
  double total_source = 0.0;
  std::array<double, N_DELAYED_GROUPS> total_delayed_source {};

  for (int g_in = 0; g_in < n_groups; ++g_in) {
    for (int g_birth = 0; g_birth < n_groups; ++g_birth) {
      const double prompt =
        data.prompt_prod_groups[prompt_source_offset(n_groups, g_in, g_birth)];
      prompt_source_by_group[g_in] += prompt;
      data.chi_prompt_groups[g_birth] += prompt;
      total_prompt_source += prompt;

      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        const double delayed = data.delayed_prod_groups[delayed_source_offset(
          n_groups, g_in, d, g_birth)];
        const size_t in_d =
          MaterialNuclearDataExtractor::delayed_offset(g_in, d);
        const size_t birth_d =
          MaterialNuclearDataExtractor::delayed_offset(g_birth, d);
        delayed_source_by_group[in_d] += delayed;
        data.chi_delayed_groups[birth_d] += delayed;
        total_delayed_source[d] += delayed;
      }
    }
  }

  for (int g = 0; g < n_groups; ++g) {
    total_source_by_group[g] = prompt_source_by_group[g];
    for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
      total_source_by_group[g] += delayed_source_by_group
        [MaterialNuclearDataExtractor::delayed_offset(g, d)];
    }
    total_source += total_source_by_group[g];

    const double sigma = data.sigma_f_groups[g];
    if (sigma <= 0.0)
      continue;

    data.nu_prompt_groups[g] = prompt_source_by_group[g] / sigma;
    data.nu_total_groups[g] = total_source_by_group[g] / sigma;
    for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
      const size_t idx = MaterialNuclearDataExtractor::delayed_offset(g, d);
      data.nu_delayed_groups[idx] = delayed_source_by_group[idx] / sigma;
    }
  }

  if (total_prompt_source > 0.0) {
    for (double& value : data.chi_prompt_groups) {
      value /= total_prompt_source;
    }
  }
  data.chi_prompt = total_prompt_source > 0.0 ? 1.0 : 0.0;

  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    if (total_delayed_source[d] <= 0.0)
      continue;
    for (int g_birth = 0; g_birth < n_groups; ++g_birth) {
      const size_t idx =
        MaterialNuclearDataExtractor::delayed_offset(g_birth, d);
      data.chi_delayed_groups[idx] /= total_delayed_source[d];
    }
    data.chi_delayed[d] = 1.0;
  }

  std::vector<double> weights(n_groups, 1.0 / static_cast<double>(n_groups));
  if (collapse_weights.size() == static_cast<size_t>(n_groups)) {
    weights = collapse_weights;
  }

  auto collapse_per_fission = [&](const std::vector<double>& values) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int g = 0; g < n_groups; ++g) {
      const double sigma = data.sigma_f_groups[g];
      if (sigma <= 0.0)
        continue;
      const double weight = weights[g] * sigma;
      numerator += weight * values[g];
      denominator += weight;
    }
    return denominator > 0.0 ? numerator / denominator : 0.0;
  };

  data.nu_total = collapse_per_fission(data.nu_total_groups);
  data.nu_prompt = collapse_per_fission(data.nu_prompt_groups);
  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    std::vector<double> delayed_slice(n_groups, 0.0);
    for (int g = 0; g < n_groups; ++g) {
      delayed_slice[g] =
        data.nu_delayed_groups[MaterialNuclearDataExtractor::delayed_offset(
          g, d)];
    }
    data.nu_delayed[d] = collapse_per_fission(delayed_slice);
  }

  if (data.sigma_f == 0.0) {
    data.sigma_f =
      std::accumulate(data.sigma_f_groups.begin(), data.sigma_f_groups.end(),
        0.0);
  }
  data.is_fissionable = total_source > 0.0 || data.sigma_f > 0.0;
}

//------------------------------------------------------------------------------
// 主入口：提取完整材料核数据
//------------------------------------------------------------------------------
MaterialNuclearData MaterialNuclearDataExtractor::extract(int material_id) const
{
  MaterialNuclearData data;
  data.material_id = material_id;

  // 查找材料
  if (model::materials.empty()) {
    warning("Material library is empty");
    return data;
  }

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
    return data;
  }

  // 1. 计算裂变截面
  if (!compute_fission_cross_sections(
        material_id, data.sigma_f_groups, data.sigma_f)) {
    return data;
  }

  // 2. Build strict source production matrices for beta_eff.
  if (!compute_fission_source_terms(
        material_id, data.prompt_prod_groups, data.delayed_prod_groups)) {
    fatal_error("Failed to build strict fission source matrices for material " +
                std::to_string(material_id) + ".");
  }

  // 3. Keep legacy diagnostic fields consistent with the strict matrices.
  derive_diagnostic_nuclear_data_from_source_terms(
    data, n_groups_, flux_weights_);

  return data;
}

//------------------------------------------------------------------------------
// 计算材料的多群裂变截面 Σ_f,g
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_fission_cross_sections(
  int material_id, std::vector<double>& sigma_f_groups,
  double& sigma_f_total) const
{
  // 初始化输出
  sigma_f_groups.assign(n_groups_, 0.0);
  sigma_f_total = 0.0;

  // 查找材料
  Material* mat = nullptr;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->id() == material_id) {
      mat = mat_ptr.get();
      break;
    }
  }
  if (!mat || !mat->fissionable()) {
    return false;
  }

  auto nuclides = mat->nuclides();
  auto densities = mat->densities();

  // 遍历所有裂变核素
  for (size_t i = 0; i < nuclides.size(); ++i) {
    int nuc_idx = nuclides[i];
    double atom_density = densities[i];

    if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
      continue;
    }

    const auto& nuc = data::nuclides[nuc_idx];
    if (!nuc || !nuc->fissionable_) {
      continue;
    }

    // 计算核素级微观截面
    std::vector<double> micro_sigma(n_groups_, 0.0);
    if (!compute_nuclide_cross_sections(
          nuc_idx, ref_temperature_, micro_sigma)) {
      continue;
    }

    // 累加宏观截面
    for (int g = 0; g < n_groups_; ++g) {
      double macro = atom_density * micro_sigma[g];
      sigma_f_groups[g] += macro;
      sigma_f_total += macro;
    }
  }

  return sigma_f_total > 0.0;
}

//------------------------------------------------------------------------------
// 计算单个核素的多群裂变截面
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_nuclide_cross_sections(
  int nuc_idx, double temperature, std::vector<double>& sigma_f_groups) const
{
  sigma_f_groups.assign(n_groups_, 0.0);

  if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
    return false;
  }

  const auto& nuc = data::nuclides[nuc_idx];
  if (!nuc || !nuc->fissionable_) {
    return false;
  }

  const int MT_FISSION = 18;
  std::vector<double> group_flux(n_groups_, 0.0);

  for (int g = 0; g < n_groups_; ++g) {
    std::fill(group_flux.begin(), group_flux.end(), 0.0);
    group_flux[g] = 1.0;

    try {
      double sigma =
        nuc->collapse_rate(MT_FISSION, temperature, energy_edges_, group_flux);
      sigma_f_groups[g] = std::max(0.0, sigma);
    } catch (...) {
      sigma_f_groups[g] = 0.0;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// 计算材料的多群中子产额
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_neutron_yields(int material_id,
  std::vector<double>& nu_total_groups, std::vector<double>& nu_prompt_groups,
  std::vector<double>& nu_delayed_groups, double& nu_total, double& nu_prompt,
  std::array<double, 8>& nu_delayed) const
{
  const size_t delayed_size = static_cast<size_t>(n_groups_) * N_DELAYED_GROUPS;

  // 初始化输出
  nu_total_groups.assign(n_groups_, 0.0);
  nu_prompt_groups.assign(n_groups_, 0.0);
  nu_delayed_groups.assign(delayed_size, 0.0);
  nu_total = 0.0;
  nu_prompt = 0.0;
  nu_delayed.fill(0.0);

  // 查找材料
  Material* mat = nullptr;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->id() == material_id) {
      mat = mat_ptr.get();
      break;
    }
  }
  if (!mat || !mat->fissionable()) {
    return false;
  }

  auto nuclides = mat->nuclides();
  auto densities = mat->densities();

  // 累加器
  std::vector<double> nu_total_accum(n_groups_, 0.0);
  std::vector<double> nu_prompt_accum(n_groups_, 0.0);
  std::vector<double> nu_delayed_accum(delayed_size, 0.0);
  std::vector<double> fission_density(n_groups_, 0.0);

  for (size_t i = 0; i < nuclides.size(); ++i) {
    int nuc_idx = nuclides[i];
    double atom_density = densities[i];

    if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
      continue;
    }

    const auto& nuc = data::nuclides[nuc_idx];
    if (!nuc || !nuc->fissionable_) {
      continue;
    }

    // 先计算核素的微观截面
    std::vector<double> micro_sigma(n_groups_, 0.0);
    compute_nuclide_cross_sections(nuc_idx, ref_temperature_, micro_sigma);

    // 计算核素的中子产额
    std::vector<double> nuc_nu_total(n_groups_, 0.0);
    std::vector<double> nuc_nu_prompt(n_groups_, 0.0);
    std::vector<std::array<double, 8>> nuc_nu_delayed(n_groups_);
    for (int g = 0; g < n_groups_; ++g) {
      nuc_nu_delayed[g].fill(0.0);
    }

    compute_nuclide_yields(nuc_idx, ref_temperature_, micro_sigma, nuc_nu_total,
      nuc_nu_prompt, nuc_nu_delayed);

    // 累加（按宏观 Σ_f 加权）
    for (int g = 0; g < n_groups_; ++g) {
      double macro_sigma = atom_density * micro_sigma[g];
      if (macro_sigma <= 0.0)
        continue;

      fission_density[g] += macro_sigma;
      nu_total_accum[g] += macro_sigma * nuc_nu_total[g];
      nu_prompt_accum[g] += macro_sigma * nuc_nu_prompt[g];

      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        nu_delayed_accum[delayed_offset(g, d)] +=
          macro_sigma * nuc_nu_delayed[g][d];
      }
    }
  }

  // 归一化：除以裂变密度
  for (int g = 0; g < n_groups_; ++g) {
    if (fission_density[g] > 0.0) {
      nu_total_groups[g] = nu_total_accum[g] / fission_density[g];
      nu_prompt_groups[g] = nu_prompt_accum[g] / fission_density[g];
      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        nu_delayed_groups[delayed_offset(g, d)] =
          nu_delayed_accum[delayed_offset(g, d)] / fission_density[g];
      }
    }
  }

  // 计算标量（通量加权）
  nu_total = flux_weighted_average(nu_total_groups, fission_density);
  nu_prompt = flux_weighted_average(nu_prompt_groups, fission_density);

  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    std::vector<double> delayed_slice(n_groups_);
    for (int g = 0; g < n_groups_; ++g) {
      delayed_slice[g] = nu_delayed_groups[delayed_offset(g, d)];
    }
    nu_delayed[d] = flux_weighted_average(delayed_slice, fission_density);
  }

  return true;
}

//------------------------------------------------------------------------------
// 计算单个核素的多群中子产额（使用 collapse_rate_weighted）
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_nuclide_yields(int nuc_idx,
  double temperature, const std::vector<double>& sigma_f_groups,
  std::vector<double>& nu_total_groups, std::vector<double>& nu_prompt_groups,
  std::vector<std::array<double, 8>>& nu_delayed_groups) const
{
  nu_total_groups.assign(n_groups_, 0.0);
  nu_prompt_groups.assign(n_groups_, 0.0);
  nu_delayed_groups.resize(n_groups_);
  for (int g = 0; g < n_groups_; ++g) {
    nu_delayed_groups[g].fill(0.0);
  }

  if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
    return false;
  }

  const auto& nuc = data::nuclides[nuc_idx];
  if (!nuc || !nuc->fissionable_) {
    return false;
  }

  const int MT_FISSION = 18;

  // 提取 delayed yield 函数
  std::array<const Function1D*, N_DELAYED_GROUPS> delayed_yield_funcs {};
  delayed_yield_funcs.fill(nullptr);
  if (!nuc->fission_rx_.empty()) {
    const auto* fission_rx = nuc->fission_rx_[0];
    int delayed_idx = 0;
    for (size_t i_prod = 1;
      i_prod < fission_rx->products_.size() && delayed_idx < N_DELAYED_GROUPS;
      ++i_prod) {
      const auto& product = fission_rx->products_[i_prod];
      if (product.particle_ == ParticleType::neutron &&
          product.emission_mode_ == ReactionProduct::EmissionMode::delayed) {
        delayed_yield_funcs[delayed_idx] = product.yield_.get();
        delayed_idx++;
      }
    }
  }

  // 定义权重函数
  std::function<double(double)> nu_total_fn = [&nuc](double E) {
    return nuc->nu(E, ReactionProduct::EmissionMode::total, 0);
  };
  std::function<double(double)> nu_prompt_fn = [&nuc](double E) {
    return nuc->nu(E, ReactionProduct::EmissionMode::prompt, 0);
  };
  std::array<std::function<double(double)>, N_DELAYED_GROUPS> delayed_fns;
  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    const Function1D* yield_fn = delayed_yield_funcs[d];
    delayed_fns[d] = [yield_fn](double E) -> double {
      if (!yield_fn)
        return 0.0;
      double val = (*yield_fn)(E);
      return val > 0.0 ? val : 0.0;
    };
  }

  // 逐群计算
  std::vector<double> group_flux(n_groups_, 0.0);

  for (int g = 0; g < n_groups_; ++g) {
    double sigma_avg = sigma_f_groups[g];
    if (sigma_avg <= 0.0)
      continue;

    std::fill(group_flux.begin(), group_flux.end(), 0.0);
    group_flux[g] = 1.0;

    try {
      double nu_total_weighted = nuc->collapse_rate_weighted(
        MT_FISSION, temperature, energy_edges_, group_flux, nu_total_fn);
      nu_total_groups[g] = nu_total_weighted / sigma_avg;

      double nu_prompt_weighted = nuc->collapse_rate_weighted(
        MT_FISSION, temperature, energy_edges_, group_flux, nu_prompt_fn);
      nu_prompt_groups[g] = nu_prompt_weighted / sigma_avg;

      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        double nu_delayed_weighted = nuc->collapse_rate_weighted(
          MT_FISSION, temperature, energy_edges_, group_flux, delayed_fns[d]);
        nu_delayed_groups[g][d] = nu_delayed_weighted / sigma_avg;
      }
    } catch (...) {
      // 失败时保持 0
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// 计算材料的多群裂变能谱
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_fission_spectra(int material_id,
  std::vector<double>& chi_prompt_groups,
  std::vector<double>& chi_delayed_groups) const
{
  const size_t delayed_size = static_cast<size_t>(n_groups_) * N_DELAYED_GROUPS;

  chi_prompt_groups.assign(n_groups_, 0.0);
  chi_delayed_groups.assign(delayed_size, 0.0);

  // 查找材料
  Material* mat = nullptr;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->id() == material_id) {
      mat = mat_ptr.get();
      break;
    }
  }
  if (!mat || !mat->fissionable()) {
    return false;
  }

  auto nuclides = mat->nuclides();
  auto densities = mat->densities();

  // 累加器
  std::vector<double> chi_prompt_acc(n_groups_, 0.0);
  std::vector<double> chi_delayed_acc(delayed_size, 0.0);
  bool prompt_sampled = false;
  std::array<bool, 8> delayed_sampled {};
  delayed_sampled.fill(false);
  double total_fission_density = 0.0;

  for (size_t i = 0; i < nuclides.size(); ++i) {
    int nuc_idx = nuclides[i];
    double atom_density = densities[i];

    if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
      continue;
    }

    const auto& nuc = data::nuclides[nuc_idx];
    if (!nuc || !nuc->fissionable_) {
      continue;
    }

    // 计算核素宏观截面（用于加权）
    std::vector<double> micro_sigma(n_groups_, 0.0);
    compute_nuclide_cross_sections(nuc_idx, ref_temperature_, micro_sigma);
    double macro_sigma_f = 0.0;
    for (int g = 0; g < n_groups_; ++g) {
      macro_sigma_f += atom_density * micro_sigma[g];
    }

    if (macro_sigma_f <= 0.0)
      continue;

    // 抽样能谱
    sample_nuclide_spectra(nuc_idx, material_id, macro_sigma_f, chi_prompt_acc,
      chi_delayed_acc, prompt_sampled, delayed_sampled);

    total_fission_density += macro_sigma_f;
  }

  // 归一化
  if (total_fission_density > 0.0) {
    for (double& val : chi_prompt_acc) {
      val /= total_fission_density;
    }
    for (double& val : chi_delayed_acc) {
      val /= total_fission_density;
    }
  }

  // 兜底：未采样的谱使用均匀分布
  double uniform = 1.0 / n_groups_;
  if (!prompt_sampled || total_fission_density <= 0.0) {
    std::fill(chi_prompt_acc.begin(), chi_prompt_acc.end(), uniform);
  }
  for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
    if (!delayed_sampled[d] || total_fission_density <= 0.0) {
      for (int g = 0; g < n_groups_; ++g) {
        chi_delayed_acc[delayed_offset(g, d)] = uniform;
      }
    }
  }

  chi_prompt_groups = chi_prompt_acc;
  chi_delayed_groups = chi_delayed_acc;

  return true;
}

//------------------------------------------------------------------------------
// 对单个核素进行能谱抽样
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::compute_fission_source_terms(int material_id,
  std::vector<double>& prompt_prod_groups,
  std::vector<double>& delayed_prod_groups) const
{
  const size_t prompt_size =
    static_cast<size_t>(n_groups_) * static_cast<size_t>(n_groups_);
  const size_t delayed_size = static_cast<size_t>(n_groups_) *
                              N_DELAYED_GROUPS *
                              static_cast<size_t>(n_groups_);
  prompt_prod_groups.assign(prompt_size, 0.0);
  delayed_prod_groups.assign(delayed_size, 0.0);

  Material* mat = nullptr;
  for (const auto& mat_ptr : model::materials) {
    if (mat_ptr && mat_ptr->id() == material_id) {
      mat = mat_ptr.get();
      break;
    }
  }
  if (!mat || !mat->fissionable()) {
    return false;
  }

  auto nuclides = mat->nuclides();
  auto densities = mat->densities();

  bool has_positive_source = false;
  bool missing_prompt_spectrum = false;
  bool missing_delayed_spectrum = false;

  for (size_t i = 0; i < nuclides.size(); ++i) {
    const int nuc_idx = nuclides[i];
    const double atom_density = densities[i];

    if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
      continue;
    }

    const auto& nuc = data::nuclides[nuc_idx];
    if (!nuc || !nuc->fissionable_) {
      continue;
    }

    std::vector<double> micro_sigma(n_groups_, 0.0);
    if (!compute_nuclide_cross_sections(
          nuc_idx, ref_temperature_, micro_sigma)) {
      warning("Missing fission cross section data for nuclide " +
              nuc->name_ + " in material " + std::to_string(material_id) +
              "; skipping this nuclide in strict source matrices.");
      continue;
    }

    std::vector<double> nuc_nu_total(n_groups_, 0.0);
    std::vector<double> nuc_nu_prompt(n_groups_, 0.0);
    std::vector<std::array<double, 8>> nuc_nu_delayed(n_groups_);
    for (int g = 0; g < n_groups_; ++g) {
      nuc_nu_delayed[g].fill(0.0);
    }
    if (!compute_nuclide_yields(nuc_idx, ref_temperature_, micro_sigma,
          nuc_nu_total, nuc_nu_prompt, nuc_nu_delayed)) {
      warning("Missing neutron yield data for nuclide " + nuc->name_ +
              " in material " + std::to_string(material_id) +
              "; skipping this nuclide in strict source matrices.");
      continue;
    }

    std::vector<double> chi_prompt_matrix;
    std::vector<double> chi_delayed_matrix;
    std::vector<bool> prompt_sampled;
    std::vector<std::array<bool, 8>> delayed_sampled;
    sample_nuclide_spectrum_matrices(nuc_idx, material_id, chi_prompt_matrix,
      chi_delayed_matrix, prompt_sampled, delayed_sampled);

    for (int g_in = 0; g_in < n_groups_; ++g_in) {
      const double macro_sigma = atom_density * micro_sigma[g_in];
      if (macro_sigma <= 0.0)
        continue;

      const double prompt_yield = nuc_nu_prompt[g_in];
      if (prompt_yield > 0.0) {
        if (g_in >= static_cast<int>(prompt_sampled.size()) ||
            !prompt_sampled[g_in]) {
          missing_prompt_spectrum = true;
        } else {
          for (int g_birth = 0; g_birth < n_groups_; ++g_birth) {
            const double chi = chi_prompt_matrix[prompt_source_offset(
              n_groups_, g_in, g_birth)];
            const double contribution = macro_sigma * prompt_yield * chi;
            prompt_prod_groups[prompt_source_offset(
              n_groups_, g_in, g_birth)] += contribution;
            has_positive_source = has_positive_source || contribution > 0.0;
          }
        }
      }

      for (int d = 0; d < N_DELAYED_GROUPS; ++d) {
        const double delayed_yield = nuc_nu_delayed[g_in][d];
        if (delayed_yield <= 0.0)
          continue;
        if (g_in >= static_cast<int>(delayed_sampled.size()) ||
            !delayed_sampled[g_in][d]) {
          missing_delayed_spectrum = true;
          continue;
        }
        for (int g_birth = 0; g_birth < n_groups_; ++g_birth) {
          const double chi = chi_delayed_matrix[delayed_source_offset(
            n_groups_, g_in, d, g_birth)];
          const double contribution = macro_sigma * delayed_yield * chi;
          delayed_prod_groups[delayed_source_offset(
            n_groups_, g_in, d, g_birth)] += contribution;
          has_positive_source = has_positive_source || contribution > 0.0;
        }
      }
    }
  }

  if (missing_prompt_spectrum) {
    warning("One or more fissionable nuclides in material " +
            std::to_string(material_id) +
            " had positive prompt yield but no sampled prompt spectrum; "
            "those prompt source terms were left zero.");
  }
  if (missing_delayed_spectrum) {
    warning("One or more fissionable nuclides in material " +
            std::to_string(material_id) +
            " had positive delayed yield but no sampled delayed spectrum; "
            "those delayed source terms were left zero.");
  }

  return has_positive_source;
}

//------------------------------------------------------------------------------

bool MaterialNuclearDataExtractor::sample_nuclide_spectra(int nuc_idx,
  int material_id, double macro_sigma_f, std::vector<double>& chi_prompt_acc,
  std::vector<double>& chi_delayed_acc, bool& prompt_sampled,
  std::array<bool, 8>& delayed_sampled) const
{
  if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
    return false;
  }

  const auto& nuc = data::nuclides[nuc_idx];
  if (!nuc || nuc->fission_rx_.empty()) {
    return false;
  }

  const auto& fission = nuc->fission_rx_[0];

  // 抽样辅助 lambda
  auto sample_spectrum = [&](const ReactionProduct& product,
                           std::vector<double>& spectrum, uint64_t seed_base) {
    if (product.distribution_.empty())
      return false;

    uint64_t seed = seed_base;
    spectrum.assign(n_groups_, 0.0);

    for (int s = 0; s < n_spectrum_samples_; ++s) {
      double E_out = 0.0, mu = 0.0;
      product.distribution_[0]->sample(ref_energy_, E_out, mu, &seed);
      int g_idx = group_index_from_energy(E_out);
      if (g_idx >= 0 && g_idx < n_groups_) {
        spectrum[g_idx] += 1.0;
      }
    }

    double total = std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
    if (total <= 0.0) {
      double uniform = 1.0 / n_groups_;
      std::fill(spectrum.begin(), spectrum.end(), uniform);
      return false;
    }
    for (double& val : spectrum) {
      val /= total;
    }
    return true;
  };

  // 采样瞬发谱
  for (const auto& product : fission->products_) {
    if (product.particle_ == ParticleType::neutron &&
        product.emission_mode_ == ReactionProduct::EmissionMode::prompt) {
      std::vector<double> prompt_spectrum;
      uint64_t seed = 0x9e3779b97f4a7c15ULL ^
                      static_cast<uint64_t>(material_id) ^
                      (static_cast<uint64_t>(nuc_idx) << 16);
      if (sample_spectrum(product, prompt_spectrum, seed)) {
        for (int g = 0; g < n_groups_; ++g) {
          chi_prompt_acc[g] += macro_sigma_f * prompt_spectrum[g];
        }
        prompt_sampled = true;
      }
      break;
    }
  }

  // 采样迟发谱
  int delayed_idx = 0;
  for (const auto& product : fission->products_) {
    if (product.particle_ == ParticleType::neutron &&
        product.emission_mode_ == ReactionProduct::EmissionMode::delayed) {
      std::vector<double> delayed_spectrum;
      uint64_t seed = 0x6a09e667f3bcc909ULL ^
                      static_cast<uint64_t>(material_id) ^
                      (static_cast<uint64_t>(nuc_idx) << 24) ^
                      static_cast<uint64_t>(delayed_idx);
      if (sample_spectrum(product, delayed_spectrum, seed)) {
        for (int g = 0; g < n_groups_; ++g) {
          chi_delayed_acc[delayed_offset(g, delayed_idx)] +=
            macro_sigma_f * delayed_spectrum[g];
        }
        delayed_sampled[delayed_idx] = true;
      }
      delayed_idx++;
      if (delayed_idx >= N_DELAYED_GROUPS)
        break;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// 通量加权平均
//------------------------------------------------------------------------------
bool MaterialNuclearDataExtractor::sample_nuclide_spectrum_matrices(
  int nuc_idx, int material_id, std::vector<double>& chi_prompt_matrix,
  std::vector<double>& chi_delayed_matrix,
  std::vector<bool>& prompt_sampled,
  std::vector<std::array<bool, 8>>& delayed_sampled) const
{
  const size_t prompt_size =
    static_cast<size_t>(n_groups_) * static_cast<size_t>(n_groups_);
  const size_t delayed_size = static_cast<size_t>(n_groups_) *
                              N_DELAYED_GROUPS *
                              static_cast<size_t>(n_groups_);
  chi_prompt_matrix.assign(prompt_size, 0.0);
  chi_delayed_matrix.assign(delayed_size, 0.0);
  prompt_sampled.assign(n_groups_, false);
  delayed_sampled.resize(n_groups_);
  for (int g = 0; g < n_groups_; ++g) {
    delayed_sampled[g].fill(false);
  }

  if (nuc_idx < 0 || nuc_idx >= static_cast<int>(data::nuclides.size())) {
    return false;
  }

  const auto& nuc = data::nuclides[nuc_idx];
  if (!nuc || nuc->fission_rx_.empty()) {
    return false;
  }

  const auto& fission = nuc->fission_rx_[0];

  auto sample_spectrum = [&](const ReactionProduct& product, int g_in,
                           std::vector<double>& spectrum, uint64_t seed_base) {
    if (product.distribution_.empty())
      return false;

    uint64_t seed = seed_base;
    spectrum.assign(n_groups_, 0.0);
    const double incident_energy = representative_incident_energy(g_in);

    for (int s = 0; s < n_spectrum_samples_; ++s) {
      double E_out = 0.0;
      double mu = 0.0;
      product.distribution_[0]->sample(incident_energy, E_out, mu, &seed);
      const int g_birth = group_index_from_energy(E_out);
      if (g_birth >= 0 && g_birth < n_groups_) {
        spectrum[g_birth] += 1.0;
      }
    }

    const double total =
      std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
    if (total <= 0.0)
      return false;

    for (double& value : spectrum) {
      value /= total;
    }
    return true;
  };

  bool any_sampled = false;

  for (const auto& product : fission->products_) {
    if (product.particle_ == ParticleType::neutron &&
        product.emission_mode_ == ReactionProduct::EmissionMode::prompt) {
      for (int g_in = 0; g_in < n_groups_; ++g_in) {
        std::vector<double> spectrum;
        const uint64_t seed = 0x9e3779b97f4a7c15ULL ^
                              static_cast<uint64_t>(material_id) ^
                              (static_cast<uint64_t>(nuc_idx) << 16) ^
                              (static_cast<uint64_t>(g_in) << 32);
        if (sample_spectrum(product, g_in, spectrum, seed)) {
          for (int g_birth = 0; g_birth < n_groups_; ++g_birth) {
            chi_prompt_matrix[prompt_source_offset(
              n_groups_, g_in, g_birth)] = spectrum[g_birth];
          }
          prompt_sampled[g_in] = true;
          any_sampled = true;
        }
      }
      break;
    }
  }

  int delayed_idx = 0;
  for (const auto& product : fission->products_) {
    if (product.particle_ == ParticleType::neutron &&
        product.emission_mode_ == ReactionProduct::EmissionMode::delayed) {
      for (int g_in = 0; g_in < n_groups_; ++g_in) {
        std::vector<double> spectrum;
        const uint64_t seed = 0x6a09e667f3bcc909ULL ^
                              static_cast<uint64_t>(material_id) ^
                              (static_cast<uint64_t>(nuc_idx) << 24) ^
                              (static_cast<uint64_t>(delayed_idx) << 40) ^
                              static_cast<uint64_t>(g_in);
        if (sample_spectrum(product, g_in, spectrum, seed)) {
          for (int g_birth = 0; g_birth < n_groups_; ++g_birth) {
            chi_delayed_matrix[delayed_source_offset(
              n_groups_, g_in, delayed_idx, g_birth)] = spectrum[g_birth];
          }
          delayed_sampled[g_in][delayed_idx] = true;
          any_sampled = true;
        }
      }
      delayed_idx++;
      if (delayed_idx >= N_DELAYED_GROUPS)
        break;
    }
  }

  return any_sampled;
}

double MaterialNuclearDataExtractor::representative_incident_energy(
  int energy_group) const
{
  if (energy_group < 0 || energy_group >= n_groups_ ||
      energy_edges_.size() != static_cast<size_t>(n_groups_ + 1)) {
    return ref_energy_;
  }

  const double e0 = energy_edges_[energy_group];
  const double e1 = energy_edges_[energy_group + 1];
  const double lo = std::min(e0, e1);
  const double hi = std::max(e0, e1);

  if (lo > 0.0 && hi > 0.0) {
    return std::sqrt(lo * hi);
  }
  const double mid = 0.5 * (lo + hi);
  return mid > 0.0 ? mid : ref_energy_;
}

//------------------------------------------------------------------------------

double MaterialNuclearDataExtractor::flux_weighted_average(
  const std::vector<double>& per_group,
  const std::vector<double>& fission_density_groups) const
{
  double accum = 0.0;
  double used_weight = 0.0;

  for (int g = 0; g < n_groups_; ++g) {
    double weight =
      (g < static_cast<int>(flux_weights_.size())) ? flux_weights_[g] : 0.0;
    if (weight <= 0.0)
      continue;
    if (g >= static_cast<int>(fission_density_groups.size()) ||
        fission_density_groups[g] <= 0.0) {
      continue;
    }
    accum += weight * per_group[g];
    used_weight += weight;
  }

  if (used_weight > 0.0) {
    return accum / used_weight;
  }
  if (!per_group.empty()) {
    return per_group.front();
  }
  return 0.0;
}

//------------------------------------------------------------------------------
// 根据能量返回所属能群索引
// 自动检测能量边界排列方向（升序或降序）
//------------------------------------------------------------------------------
int MaterialNuclearDataExtractor::group_index_from_energy(
  double energy_eV) const
{
  if (energy_edges_.size() < 2)
    return 0;

  // 检测能量边界排列方向
  bool ascending = (energy_edges_.front() < energy_edges_.back());

  if (ascending) {
    // 能量边界升序排列（低能在前）: edges[g] <= E < edges[g+1]
    if (energy_eV <= energy_edges_.front())
      return 0;
    if (energy_eV >= energy_edges_.back())
      return n_groups_ - 1;

    for (int g = 0; g < n_groups_; ++g) {
      if (energy_eV >= energy_edges_[g] && energy_eV < energy_edges_[g + 1]) {
        return g;
      }
    }
  } else {
    // 能量边界降序排列（高能在前）: edges[g+1] < E <= edges[g]
    if (energy_eV >= energy_edges_.front())
      return 0;
    if (energy_eV <= energy_edges_.back())
      return n_groups_ - 1;

    for (int g = 0; g < n_groups_; ++g) {
      if (energy_eV <= energy_edges_[g] && energy_eV > energy_edges_[g + 1]) {
        return g;
      }
    }
  }
  return n_groups_ - 1;
}

//------------------------------------------------------------------------------
// 获取平铺数组中的延迟群索引
//------------------------------------------------------------------------------
size_t MaterialNuclearDataExtractor::delayed_offset(
  int energy_group, int delayed_group)
{
  return static_cast<size_t>(energy_group) * N_DELAYED_GROUPS + delayed_group;
}

} // namespace openmc
