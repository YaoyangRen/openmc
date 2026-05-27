#ifndef OPENMC_BETA_EFFECTIVE_H
#define OPENMC_BETA_EFFECTIVE_H

#include "openmc/material_nuclear_data.h"

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

//==============================================================================
//! 单种方法的 β_eff 结果集
//==============================================================================
struct MethodResult {
  std::array<double, N_DELAYED_GROUPS> beta_i {};     //!< 各组 β_i,eff
  std::array<double, N_DELAYED_GROUPS> numerators {}; //!< 各组分子
  double denominator {0.0};                           //!< 分母
  double beta_total {0.0};                            //!< 总 β_eff
  bool available {false};                             //!< 方法是否成功计算
};

// Forward declarations
struct Position;

//==============================================================================
//! 多点采样加权模式枚举 (Phase 2.3)
//==============================================================================
enum class WeightingMode {
  VOLUME_WEIGHTED,       //!< 体积分数加权
  REACTION_RATE_WEIGHTED //!< 反应率加权 (物理更准确)
};

//==============================================================================
//! 有效缓发中子份额计算类
//!
//! Current main path uses fission neutron birth-energy source importance:
//!
//! D = sum_c dV sum_gin sum_gb phi_gin Sigma_f,gin
//!     [nu_p,gin chi_p,gb I*(c,gb)
//!      + sum_k nu_d,k,gin chi_d,k,gb I*(c,gb)]
//! W_t(c,gin) is the bracketed birth-spectrum-folded total source importance.
//! N_k = sum_c dV sum_gin phi_gin Sigma_f,gin
//!       (nu_d,k,gin / nu_t,gin) W_t(c,gin)
//!
//! I*(c,gb) is read from fission_matrix.h5/adjoint_source_grouped. The
//! response-collision-energy Method E path is retained as method comparison.
//==============================================================================
class BetaEffective {
public:
  //! 构造函数
  //! \param n_sample_points 每个网格单元的采样点数 (1, 8, 或 27)
  //! \param weighting_mode 多材料加权模式
  //! \param ref_energy 参考中子能量 (eV)
  //! \param ref_temperature 参考温度 (K)
  explicit BetaEffective(int n_sample_points = 27,
    WeightingMode weighting_mode = WeightingMode::REACTION_RATE_WEIGHTED,
    double ref_energy = 2.0e6, double ref_temperature = 293.6)
    : n_sample_points_(n_sample_points), weighting_mode_(weighting_mode),
      ref_energy_(ref_energy), ref_temperature_(ref_temperature)
  {}

  //! 从 HDF5 文件计算 β_eff
  //! \param flux_file 正向通量文件 (flux_mesh.h5)
  //! \param adjoint_flux_file response-weighted importance 文件 (adjoint_flux.h5)
  //! \param output_file 输出文件 (beta_eff.h5)
  void compute_from_files(const std::string& flux_file,
    const std::string& adjoint_flux_file, const std::string& output_file,
    const std::string& fission_matrix_file = "fission_matrix.h5");

  //! 获取特定缓发群的 β_i,eff
  //! \param group 缓发群索引 (0-7)
  double get_beta_i(int group) const;

  //! 获取总 β_eff
  double get_beta_total() const { return beta_total_; }

  //! 获取所有缓发群的 β_i,eff
  const std::array<double, N_DELAYED_GROUPS>& get_all_beta_i() const
  {
    return beta_i_;
  }

private:
  //==========================================================================
  // Birth-energy source-state importance method (main result)
  //==========================================================================

  //! Main denominator: total birth source folded with I*(cell,g_birth)
  double compute_denominator_birth_spectrum(
    const std::unordered_map<int, double>& flux, double volume);

  //! Main numerator: total birth-source importance apportioned by nu_d,k/nu_t
  double compute_delayed_numerator_birth_spectrum(int group,
    const std::unordered_map<int, double>& flux, double volume);

  //==========================================================================
  // Method E: upstream family-resolved response-weighted importance
  //==========================================================================

  //! 方法 E 分母: D = Σ_cell ΔV × Σ_g I_total(cell,g) × F_total(cell,g)
  //! I_total = I_prompt + Σ_k I_delayed_k (来自上游族解析 importance)
  double compute_denominator_upstream_family(
    const std::unordered_map<int, double>& flux, double volume);

  //! 方法 E 分子: N_k = Σ_cell ΔV × Σ_g I_delayed_k(cell,g) × F_total(cell,g)
  //! I_delayed_k 来自上游族解析 importance (传递函数 family 轴)
  double compute_delayed_numerator_upstream_family(int group,
    const std::unordered_map<int, double>& flux, double volume);

  //! 输出详细诊断信息
  void print_diagnostic_info(const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 从材料库提取核数据
  //! \param material_id OpenMC 材料 ID
  MaterialNuclearData extract_material_nuclear_data(int material_id) const;

  //! 构建单元-材料映射
  //! 通过几何查询确定每个网格单元的实际材料
  //! 把稀疏网格上的每个非零单元映射到其真实几何材料组成，
  //! 并预先提取/混合好材料核数据（Σ_f、ν_total/ν_delayed、χ 等）。
  //! 这样计算分子/分母时就能直接查 cell_nuclear_data_，避免在 β_eff
  //! 主循环里重复几何查找或核数据提取。
  void build_cell_material_map(const std::unordered_map<int, double>& flux);

  //! 生成规则网格采样位置
  //! \param lower 单元左下角坐标
  //! \param upper 单元右上角坐标
  //! \param point_idx 采样点索引 (0-7 for 8点, 0-26 for 27点)
  //! \param n_points 总采样点数 (8 或 27)
  Position sample_cell_point(const Position& lower, const Position& upper,
    int point_idx, int n_points) const;

  //! 计算多材料单元的加权核数据
  //! \param material_counts 材料 -> 命中次数的映射
  //! \param total_samples 总采样点数
  //! \param material_cache 材料核数据缓存（线程安全）
  MaterialNuclearData compute_weighted_nuclear_data(
    const std::unordered_map<int, int>& material_counts, int total_samples,
    const std::unordered_map<int, MaterialNuclearData>& material_cache) const;

  //! Kahan 求和算法（提高数值精度）
  //! \param values 待求和的数值向量
  //! \return 高精度求和结果
  double kahan_sum(const std::vector<double>& values) const;

  //! 从网格索引计算空间位置 (旧方法，保留兼容性)
  std::array<double, 3> mesh_index_to_position(int mesh_idx) const;

  //! 从一维网格索引计算三维索引
  //! \param cell_idx 一维单元索引 (0-based)
  //! \return 三维网格索引 {i, j, k}
  std::array<int, 3> get_grid_indices(int cell_idx) const;

  //! 从网格单元索引计算中心位置
  //! \param cell_idx 一维单元索引
  //! \return 单元中心的全局坐标 (用于几何查询)
  Position grid_index_to_position(int cell_idx) const;

  //! 从 HDF5 读取稀疏通量数据
  void read_flux_data(const std::string& filename,
    std::unordered_map<int, double>& flux_map, std::array<int, 3>& shape,
    double& pitch);

  //! 从 HDF5 读取稀疏 response-weighted importance 数据
  void read_adjoint_flux_data(const std::string& filename,
    std::unordered_map<int, double>& adjoint_flux_map,
    std::array<int, 3>& shape, double& pitch);

  //! Read source-state importance I*(cell,g_birth) from fission_matrix.h5
  void read_fission_adjoint_source_data(const std::string& filename,
    const std::array<int, 3>& expected_shape, double expected_pitch);

  //! 校验 (并必要时设置) 通量/importance 的能群一致性
  void validate_group_metadata();

  //! 根据正向通量文件构建能谱折算权重
  //! 对所有cell累加各能群通量，归一化后作为后续求解发射谱的权重
  void initialize_flux_spectrum_weights();

  //! 根据当前能量网格返回能群索引 (落在区间外时夹紧)
  int group_index_from_energy(double energy_eV) const;

  //! 获取平铺数组中的延迟群索引
  size_t delayed_offset(int energy_group, int delayed_group) const;

  //! 写入结果到 HDF5 文件
  void write_to_file(const std::string& filename) const;

  // ========== 数据成员 ==========
  int n_sample_points_;          //!< 每个网格单元的采样点数
  WeightingMode weighting_mode_; //!< 多材料加权模式
  double ref_energy_;            //!< 参考中子能量 (eV)
  double ref_temperature_;       //!< 参考温度 (K)

  // 使用 material_nuclear_data.h 中定义的 N_DELAYED_GROUPS = 8

  std::array<double, N_DELAYED_GROUPS> beta_i_;     //!< 各组 β_i,eff
  double beta_total_;                               //!< 总 β_eff
  std::array<double, N_DELAYED_GROUPS> numerators_; //!< 各组分子
  double denominator_;                              //!< 分母

  // 多群通量数据缓存
  int flux_n_groups_ {1};
  int adjoint_n_groups_ {1};
  int n_energy_groups_ {1};
  std::vector<double> flux_energy_edges_;
  std::vector<double> adjoint_energy_edges_;
  std::vector<double> energy_edges_common_;
  std::unordered_map<int, std::vector<double>> flux_group_map_;
  std::unordered_map<int, std::vector<double>> adjoint_group_map_;
  bool flux_has_group_data_ {false};
  bool adjoint_has_group_data_ {false};
  std::vector<double> flux_collapse_weights_; //!< 归一化通量权重(按能群)

  // Phase 2.3: 多点采样统计
  int n_heterogeneous_cells_ {0}; //!< 包含多材料的单元数

  //! 方法 E 结果缓存
  MethodResult result_e_; //!< 方法 E: 上游族解析

  // Method E: 上游族解析 importance 缓存
  MethodResult result_birth_spectrum_; //!< Main: birth-energy source method

  bool has_upstream_family_data_ {false};
  bool upstream_family_group_data_used_ {false};
  std::unordered_map<int, double> cell_upstream_prompt_importance_;
  std::array<std::unordered_map<int, double>, N_DELAYED_GROUPS>
    cell_upstream_delayed_importance_;
  std::unordered_map<int, std::vector<double>>
    cell_upstream_prompt_importance_group_;
  std::array<std::unordered_map<int, std::vector<double>>, N_DELAYED_GROUPS>
    cell_upstream_delayed_importance_group_;

  // Method E energy-group diagnostics
  std::array<std::vector<double>, N_DELAYED_GROUPS>
    numerator_by_group_energy_;
  std::vector<double> denominator_by_group_energy_;

  // Birth-energy source-state importance I*(cell,g_birth)
  bool has_birth_adjoint_source_ {false};
  int birth_source_n_groups_ {1};
  std::vector<double> birth_source_energy_edges_;
  std::vector<double> birth_adjoint_source_grouped_;
  std::string fission_matrix_file_ {"fission_matrix.h5"};

  // 网格信息(用于验证一致性)
  std::array<int, 3> grid_shape_;
  double grid_pitch_;
  std::array<double, 3> grid_lower_left_; //!< 网格左下角坐标

  // Phase 2.2: 材料相关核数据 (几何查询)
  std::unordered_map<int, MaterialNuclearData>
    cell_nuclear_data_;                               //!< 单元核数据映射
  std::unordered_map<int, int> cell_to_material_;     //!< 单元到材料ID的映射
  std::vector<MaterialNuclearData> unique_materials_; //!< 唯一材料列表(诊断用)
};

} // namespace openmc

#endif // OPENMC_BETA_EFFECTIVE_H
