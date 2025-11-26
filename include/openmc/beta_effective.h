#ifndef OPENMC_BETA_EFFECTIVE_H
#define OPENMC_BETA_EFFECTIVE_H

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

// Forward declarations
struct Position;

//==============================================================================
//! 材料核数据结构 (Phase 2 -> Phase 3)
//! 存储单个材料的裂变相关核参数
//==============================================================================
struct MaterialNuclearData {
  double nu_total {0.0};               //!< 总裂变中子产额
  double nu_prompt {0.0};              //!< 瞬发中子产额
  std::array<double, 8> nu_delayed {}; //!< 8组缓发中子产额 (Phase 3)
  double sigma_f {0.0};                //!< 宏观裂变截面 (cm^-1)
  double chi_prompt {1.0};             //!< 瞬发中子能谱 (归一化)
  std::array<double, 8> chi_delayed {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; //!< 缓发中子能谱 (Phase 3)
  int material_id {-1};                      //!< 材料ID
  std::string material_name;                 //!< 材料名称
  bool is_fissionable {false};               //!< 是否可裂变
};

//==============================================================================
//! β_eff 计算模式枚举
//==============================================================================
enum class BetaEffMode {
  FIXED_U235,        //!< Phase 1: 固定 U-235 核数据
  MATERIAL_DEPENDENT //!< Phase 2: 材料相关核数据
};

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
  //! \param mode 计算模式 (FIXED_U235 或 MATERIAL_DEPENDENT)
  //! \param n_sample_points 每个网格单元的采样点数 (Phase 2.3: 1, 8, 或 27)
  //! \param weighting_mode 多材料加权模式 (Phase 2.3)
  //! \param ref_energy 参考中子能量 (eV, Phase 3)
  //! \param ref_temperature 参考温度 (K, Phase 3)
  explicit BetaEffective(BetaEffMode mode = BetaEffMode::MATERIAL_DEPENDENT,
    int n_sample_points = 27,
    WeightingMode weighting_mode = WeightingMode::REACTION_RATE_WEIGHTED,
    double ref_energy = 0.0253, double ref_temperature = 293.6)
    : mode_(mode), n_sample_points_(n_sample_points),
      weighting_mode_(weighting_mode), ref_energy_(ref_energy),
      ref_temperature_(ref_temperature)
  {}

  //! 从 HDF5 文件计算 β_eff
  //! \param flux_file 正向通量文件 (flux_mesh.h5)
  //! \param adjoint_flux_file 共轭通量文件 (adjoint_flux.h5)
  //! \param output_file 输出文件 (beta_eff.h5)
  void compute_from_files(const std::string& flux_file,
    const std::string& adjoint_flux_file, const std::string& output_file);

  //! 获取特定缓发群的 β_i,eff
  //! \param group 缓发群索引 (0-7, Phase 3)
  double get_beta_i(int group) const;

  //! 获取总 β_eff
  double get_beta_total() const { return beta_total_; }

  //! 获取所有缓发群的 β_i,eff
  const std::array<double, 8>& get_all_beta_i() const { return beta_i_; }

private:
  //! 计算缓发中子贡献的分子 (Phase 1)
  //! \param group 缓发群索引 (0-7, Phase 3)
  //! \param flux 正向通量 (稀疏格式)
  //! \param adjoint_flux 共轭通量 (稀疏格式)
  //! \param volume 网格单元体积
  double compute_delayed_numerator(int group,
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算总中子贡献的分母 (Phase 1)
  double compute_denominator(const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算缓发中子贡献的分子 (Phase 2: 材料相关)
  double compute_delayed_numerator_material(int group,
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算总中子贡献的分母 (Phase 2: 材料相关)
  double compute_denominator_material(
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 从材料库提取核数据 (Phase 2)
  //! \param material_id OpenMC 材料 ID
  MaterialNuclearData extract_material_nuclear_data(int material_id) const;

  //! 构建单元-材料映射 (Phase 2.2/2.3)
  //! 通过几何查询确定每个网格单元的实际材料
  void build_cell_material_map(const std::unordered_map<int, double>& flux);

  //! 生成规则网格采样位置 (Phase 2.3)
  //! \param lower 单元左下角坐标
  //! \param upper 单元右上角坐标
  //! \param point_idx 采样点索引 (0-7 for 8点, 0-26 for 27点)
  //! \param n_points 总采样点数 (8 或 27)
  Position sample_cell_point(const Position& lower, const Position& upper,
    int point_idx, int n_points) const;

  //! 计算多材料单元的加权核数据 (Phase 2.3)
  //! \param material_counts 材料 -> 命中次数的映射
  //! \param total_samples 总采样点数
  MaterialNuclearData compute_weighted_nuclear_data(
    const std::unordered_map<int, int>& material_counts,
    int total_samples) const;

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
    std::array<int, 3>& shape, double& pitch) const;

  //! 写入结果到 HDF5 文件
  void write_to_file(const std::string& filename) const;

  // ========== 数据成员 ==========
  BetaEffMode mode_;             //!< 计算模式
  int n_sample_points_;          //!< 每个网格单元的采样点数 (Phase 2.3)
  WeightingMode weighting_mode_; //!< 多材料加权模式 (Phase 2.3)
  double ref_energy_;            //!< 参考中子能量 (eV, Phase 3)
  double ref_temperature_;       //!< 参考温度 (K, Phase 3)

  std::array<double, 8> beta_i_;     //!< 各组 β_i,eff (Phase 3: 8组)
  double beta_total_;                //!< 总 β_eff
  std::array<double, 8> numerators_; //!< 各组分子(诊断用, Phase 3: 8组)
  double denominator_;               //!< 分母(诊断用)

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

  // Phase 1: 使用 U-235 热中子裂变的固定核数据
  static constexpr double NU_TOTAL = 2.43;  //!< 总中子产额
  static constexpr double NU_PROMPT = 2.42; //!< 瞬发中子产额
  static constexpr double CHI_PROMPT = 1.0; //!< 瞬发中子谱(归一化)
  static constexpr double SIGMA_F = 1.0;    //!< 相对裂变截面(归一化)

  //! U-235 的 6 组缓发中子参数
  static constexpr std::array<double, 6> NU_DELAYED = {
    0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273};

  //! 各组缓发中子谱(简化假设为相同)
  static constexpr std::array<double, 6> CHI_DELAYED = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};

} // namespace openmc

#endif // OPENMC_BETA_EFFECTIVE_H
