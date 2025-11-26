#ifndef OPENMC_BETA_EFFECTIVE_H
#define OPENMC_BETA_EFFECTIVE_H

#include <array>
#include <string>
#include <unordered_map>

namespace openmc {

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
  //! 从 HDF5 文件计算 β_eff
  //! \param flux_file 正向通量文件 (flux_mesh.h5)
  //! \param adjoint_flux_file 共轭通量文件 (adjoint_flux.h5)
  //! \param output_file 输出文件 (beta_eff.h5)
  void compute_from_files(const std::string& flux_file,
    const std::string& adjoint_flux_file, const std::string& output_file);

  //! 获取特定缓发群的 β_i,eff
  //! \param group 缓发群索引 (0-5)
  double get_beta_i(int group) const;

  //! 获取总 β_eff
  double get_beta_total() const { return beta_total_; }

  //! 获取所有缓发群的 β_i,eff
  const std::array<double, 6>& get_all_beta_i() const { return beta_i_; }

private:
  //! 计算缓发中子贡献的分子
  //! \param group 缓发群索引 (0-5)
  //! \param flux 正向通量 (稀疏格式)
  //! \param adjoint_flux 共轭通量 (稀疏格式)
  //! \param volume 网格单元体积
  double compute_delayed_numerator(int group,
    const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 计算总中子贡献的分母
  double compute_denominator(const std::unordered_map<int, double>& flux,
    const std::unordered_map<int, double>& adjoint_flux, double volume) const;

  //! 从 HDF5 读取稀疏通量数据
  void read_flux_data(const std::string& filename,
    std::unordered_map<int, double>& flux_map, std::array<int, 3>& shape,
    double& pitch) const;

  //! 从 HDF5 读取稀疏共轭通量数据
  void read_adjoint_flux_data(const std::string& filename,
    std::unordered_map<int, double>& adjoint_flux_map,
    std::array<int, 3>& shape, double& pitch) const;

  //! 写入结果到 HDF5 文件
  void write_to_file(const std::string& filename) const;

  // 数据成员
  std::array<double, 6> beta_i_;     //!< 各组 β_i,eff
  double beta_total_;                //!< 总 β_eff
  std::array<double, 6> numerators_; //!< 各组分子(诊断用)
  double denominator_;               //!< 分母(诊断用)

  // 网格信息(用于验证一致性)
  std::array<int, 3> grid_shape_;
  double grid_pitch_;

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
