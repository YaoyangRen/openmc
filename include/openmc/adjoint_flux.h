#ifndef OPENMC_ADJOINT_FLUX_H
#define OPENMC_ADJOINT_FLUX_H

#include "openmc/array.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <string>
#include <unordered_map>

namespace openmc {

//! \class AdjointFlux
//! \brief 共轭通量（伴随通量）计算器
//!
//! 物理意义：
//!   Φ†(r) = Σ_i T(i -> r) × S†(i)
//!
//!   其中：
//!   - T(i -> r): 传递函数，从源单元 i 到响应位置 r 的传递
//!   - S†(i): 伴随源（重要性函数），源单元 i 的重要性
//!   - Φ†(r): 共轭通量，位置 r 处的伴随通量
//!
//! 用途：
//!   - 扰动理论分析
//!   - 灵敏度系数计算
//!   - 重要性采样权重
//!   - 不确定性量化
//!
class AdjointFlux {
public:
  // 默认构造函数
  AdjointFlux() = default;

  // 从 HDF5 文件计算共轭通量
  // transfer_function_file: 传递函数数据文件 (transfer_function_data.h5)
  // fission_matrix_file: 裂变矩阵文件，包含伴随源 (fission_matrix.h5)
  // output_file: 输出的共轭通量文件 (adjoint_flux.h5)
  void compute_from_files(const std::string& transfer_function_file,
    const std::string& fission_matrix_file,
    const std::string& output_file = "adjoint_flux.h5");

  // 从内存中的数据计算共轭通量
  // transfer_functions: 稀疏传递函数 map<source_state, map<j_response, T>>
  //   source_state = source_cell * n_source_groups + g_source
  //   若 n_source_groups==1，source_state == source_cell（向后兼容）
  // adjoint_source: 伴随源分布向量
  //   若 n_source_groups==1：size = n_cells（标量 I*(cell)）
  //   若 n_source_groups>1 ：size = n_cells*n_source_groups（分群 I*(cell,g)）
  // n_source_groups: 源能群数（默认 1，等于标量模式）
  void compute_from_memory(
    const std::unordered_map<int, std::unordered_map<int, vector<double>>>&
      transfer_functions,
    const vector<double>& adjoint_source, const std::array<int, 3>& shape,
    const std::array<double, 3>& origin, double pitch, int n_groups = 1,
    int n_families = 1, vector<double> energy_edges = {},
    int n_source_groups = 1);

  // 获取计算结果（稀疏格式）
  const std::unordered_map<int, vector<double>>& get_adjoint_flux_sparse() const
  {
    return adjoint_flux_sparse_;
  }

  // 获取计算结果（稠密格式）
  vector<double> get_adjoint_flux_dense() const;
  vector<double> get_adjoint_flux_group_dense() const;

  // 写入 HDF5 文件
  void write_to_file(const std::string& filename = "adjoint_flux.h5");

  // 获取网格信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }
  int n_groups() const { return n_groups_; }
  const vector<double>& energy_edges() const { return energy_edges_; }

  // 获取统计信息
  size_t get_nonzero_count() const { return adjoint_flux_sparse_.size(); }
  double get_max_value() const;
  double get_total_flux() const;

  // Family-resolved accessors
  bool has_family_data() const { return has_family_data_; }
  int n_families() const { return n_families_; }
  const std::vector<std::unordered_map<int, vector<double>>>&
  get_family_adjoint_flux() const
  {
    return family_adjoint_flux_;
  }

private:
  // 共轭通量数据（稀疏存储）
  // Key: 单元索引 j, Value: Φ†(j)
  std::unordered_map<int, vector<double>> adjoint_flux_sparse_;

  // 网格参数
  std::array<int, 3> shape_ {0, 0, 0};
  std::array<double, 3> origin_ {0.0, 0.0, 0.0};
  double pitch_ {0.0};
  size_t n_cells_ {0};
  int n_groups_ {1};
  int n_families_ {1};
  bool has_family_data_ {false};
  vector<double> energy_edges_;
  vector<double> group_total_flux_;

  // Family-resolved adjoint flux: [family_index] -> (cell -> group_vector)
  std::vector<std::unordered_map<int, vector<double>>> family_adjoint_flux_;

  // 统计信息
  double max_flux_ {0.0};
  double total_flux_ {0.0};
  size_t nonzero_cells_ {0};

  std::vector<double> make_zero_group_vector() const;
};

} // namespace openmc

#endif // OPENMC_ADJOINT_FLUX_H
