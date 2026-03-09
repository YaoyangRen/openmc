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

  // 2. 计算中子产额
  if (!compute_neutron_yields(material_id, data.nu_total_groups,
        data.nu_prompt_groups, data.nu_delayed_groups, data.nu_total,
        data.nu_prompt, data.nu_delayed)) {
    return data;
  }

  // 3. 计算裂变能谱
  if (!compute_fission_spectra(
        material_id, data.chi_prompt_groups, data.chi_delayed_groups)) {
    // 能谱计算失败时使用均匀分布
    double uniform = 1.0 / n_groups_;
    data.chi_prompt_groups.assign(n_groups_, uniform);
    data.chi_delayed_groups.assign(
      static_cast<size_t>(n_groups_) * N_DELAYED_GROUPS, uniform);
  }

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
