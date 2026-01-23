#ifndef OPENMC_MATERIAL_NUCLEAR_DATA_H
#define OPENMC_MATERIAL_NUCLEAR_DATA_H

#include <array>
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
//! 常量定义
//==============================================================================
constexpr int N_DELAYED_GROUPS = 8; //!< 缓发中子群数（与 ENDF 一致）

//==============================================================================
//! 材料核数据结构
//! 存储单个材料的裂变相关核参数（支持单群和多群）
//==============================================================================
struct MaterialNuclearData {
  // --- 单能群参考值（标量，用于诊断/输出）---
  double nu_total {0.0};                //!< 总裂变中子产额
  double nu_prompt {0.0};               //!< 瞬发中子产额
  std::array<double, 8> nu_delayed {};  //!< 8组缓发中子产额
  double sigma_f {0.0};                 //!< 宏观裂变截面 (cm^-1)
  double chi_prompt {1.0};              //!< 单群瞬发能谱权重
  std::array<double, 8> chi_delayed {}; //!< 单群缓发能谱权重

  // --- 多能群扩展数据（与 flux/adjoint 能群对齐）---
  std::vector<double> sigma_f_groups;     //!< Σ_f,g (cm^-1)
  std::vector<double> nu_total_groups;    //!< ν_total,g
  std::vector<double> nu_prompt_groups;   //!< ν_prompt,g
  std::vector<double> nu_delayed_groups;  //!< ν_delayed,g,i (flat: g-major)
  std::vector<double> chi_prompt_groups;  //!< χ_prompt,g
  std::vector<double> chi_delayed_groups; //!< χ_delayed,g,i (flat)

  int material_id {-1};        //!< 材料ID
  std::string material_name;   //!< 材料名称
  bool is_fissionable {false}; //!< 是否可裂变
};

//==============================================================================
//! 核素级核数据结构（用于中间计算）
//==============================================================================
struct NuclideNuclearData {
  std::vector<double> sigma_f_groups;                   //!< 微观 σ_f,g (barn)
  std::vector<double> nu_total_groups;                  //!< ν_total,g
  std::vector<double> nu_prompt_groups;                 //!< ν_prompt,g
  std::vector<std::array<double, 8>> nu_delayed_groups; //!< ν_delayed,g,d
  std::vector<double> chi_prompt_groups;  //!< χ_prompt,g（归一化后）
  std::vector<double> chi_delayed_groups; //!< χ_delayed,g,d（flat，归一化后）
  double atom_density {0.0};              //!< 原子密度 (atom/b-cm)
  int nuclide_index {-1};                 //!< 核素索引
  std::string nuclide_name;               //!< 核素名称
  bool is_fissionable {false};
};

//==============================================================================
//! 材料核数据提取器类
//! 将材料的裂变截面、中子产额、能谱分离为独立函数
//==============================================================================
class MaterialNuclearDataExtractor {
public:
  //! 构造函数
  //! \param n_energy_groups 能群数
  //! \param energy_edges 能量边界 (eV)，长度 = n_groups + 1
  //! \param flux_weights 通量谱折合权重（归一化），用于标量 ν 的计算
  //! \param ref_energy 参考入射能量 (eV)，用于能谱抽样
  //! \param ref_temperature 参考温度 (K)
  MaterialNuclearDataExtractor(int n_energy_groups,
    const std::vector<double>& energy_edges,
    const std::vector<double>& flux_weights = {}, double ref_energy = 2.0e6,
    double ref_temperature = 293.6);

  //! 提取完整材料核数据（调用下面三个子函数）
  //! \param material_id 材料 ID
  //! \return 填充好的 MaterialNuclearData
  MaterialNuclearData extract(int material_id) const;

  //! 计算材料的多群裂变截面 Σ_f,g
  //! \param material_id 材料 ID
  //! \param[out] sigma_f_groups 各群宏观裂变截面 (cm^-1)
  //! \param[out] sigma_f_total 总宏观裂变截面
  //! \return 是否成功（非裂变材料返回 false）
  bool compute_fission_cross_sections(int material_id,
    std::vector<double>& sigma_f_groups, double& sigma_f_total) const;

  //! 计算材料的多群中子产额 ν_total,g / ν_prompt,g / ν_delayed,g,d
  //! \param material_id 材料 ID
  //! \param[out] nu_total_groups 各群总产额
  //! \param[out] nu_prompt_groups 各群瞬发产额
  //! \param[out] nu_delayed_groups 各群各迟发组产额 (flat: g*8+d)
  //! \param[out] nu_total 标量总产额（通量加权）
  //! \param[out] nu_prompt 标量瞬发产额
  //! \param[out] nu_delayed 标量各迟发组产额
  //! \return 是否成功
  bool compute_neutron_yields(int material_id,
    std::vector<double>& nu_total_groups, std::vector<double>& nu_prompt_groups,
    std::vector<double>& nu_delayed_groups, double& nu_total, double& nu_prompt,
    std::array<double, 8>& nu_delayed) const;

  //! 计算材料的多群裂变能谱 χ_prompt,g / χ_delayed,g,d
  //! \param material_id 材料 ID
  //! \param[out] chi_prompt_groups 瞬发能谱 (归一化)
  //! \param[out] chi_delayed_groups 迟发能谱 (flat: g*8+d，归一化)
  //! \return 是否成功
  bool compute_fission_spectra(int material_id,
    std::vector<double>& chi_prompt_groups,
    std::vector<double>& chi_delayed_groups) const;

  //! 根据能量返回所属能群索引
  int group_index_from_energy(double energy_eV) const;

  //! 获取平铺数组中的延迟群索引
  static size_t delayed_offset(int energy_group, int delayed_group);

private:
  //! 计算单个核素的多群裂变截面
  bool compute_nuclide_cross_sections(
    int nuc_idx, double temperature, std::vector<double>& sigma_f_groups) const;

  //! 计算单个核素的多群中子产额（σ_f 加权）
  bool compute_nuclide_yields(int nuc_idx, double temperature,
    const std::vector<double>& sigma_f_groups,
    std::vector<double>& nu_total_groups, std::vector<double>& nu_prompt_groups,
    std::vector<std::array<double, 8>>& nu_delayed_groups) const;

  //! 对单个核素的裂变产物进行能谱抽样
  bool sample_nuclide_spectra(int nuc_idx, int material_id,
    double macro_sigma_f, std::vector<double>& chi_prompt_acc,
    std::vector<double>& chi_delayed_acc, bool& prompt_sampled,
    std::array<bool, 8>& delayed_sampled) const;

  //! 通量加权平均（将群值折合为标量）
  double flux_weighted_average(const std::vector<double>& per_group,
    const std::vector<double>& fission_density_groups) const;

  // 配置参数
  int n_groups_;
  std::vector<double> energy_edges_;
  std::vector<double> flux_weights_;
  double ref_energy_;
  double ref_temperature_;
  int n_spectrum_samples_ {5000}; //!< 能谱抽样次数
};

} // namespace openmc

#endif // OPENMC_MATERIAL_NUCLEAR_DATA_H
