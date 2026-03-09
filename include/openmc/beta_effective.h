#ifndef OPENMC_BETA_EFFECTIVE_H
#define OPENMC_BETA_EFFECTIVE_H

#include "openmc/material_nuclear_data.h"

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

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
//! 计算 β_eff = Σ_i β_i,eff，其中每个缓发群的有效份额为:
//! β_i,eff = (分子_i) / (分母)
//!
//! 分子_i = ∫ φ*(r) χ_d,i ν_d,i(r) σ_f(r) φ(r) dr
//! 分母 = ∫ φ*(r) χ_p ν(r) σ_f(r) φ(r) dr
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
  //! \param adjoint_flux_file 共轭通量文件 (adjoint_flux.h5)
  //! \param output_file 输出文件 (beta_eff.h5)
  void compute_from_files(const std::string& flux_file,
    const std::string& adjoint_flux_file, const std::string& output_file);

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
  //! 计算缓发中子贡献的分子
  double compute_delayed_numerator_material(int group,
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算总中子贡献的分母 (含 chi_prompt)
  double compute_denominator_material(
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算不含 chi 的分母 (用于诊断对比)
  double compute_denominator_without_chi(
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

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

  //! 从 HDF5 读取稀疏共轭通量数据
  void read_adjoint_flux_data(const std::string& filename,
    std::unordered_map<int, double>& adjoint_flux_map,
    std::array<int, 3>& shape, double& pitch);

  //! 校验 (并必要时设置) 通量/共轭通量的能群一致性
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
