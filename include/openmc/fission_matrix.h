#ifndef OPENMC_FISSION_MATRIX_H
#define OPENMC_FISSION_MATRIX_H

#include "hdf5.h"
#include "openmc/array.h"
#include "openmc/mesh_init.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class FissionMatrix {
public:
  // 构造函数
  // energy_edges: 源中子出生能群边界（与 settings::kinetics_energy_edges 共享）
  // 若为空则退化为标量模式（n_source_groups=1），保持向后兼容
  explicit FissionMatrix(std::shared_ptr<SharedMeshGrid> grid, int max_batches,
    std::vector<double> energy_edges = {});

  // 记录源中子的产生位置
  // energy: CE 模式下源粒子出生能量 (eV)；MG 模式下置 -1 并传入 mg_group
  // mg_group: MG 模式下出生群号；CE 模式下置 -1
  void record_source_birth(const Position& r, int64_t source_particle_id,
    double energy = -1.0, int mg_group = -1);

  // 记录裂变事件（裂变位置和产生的中子数）
  void record_fission_event(
    const Position& r, double nu_fission, int64_t source_particle_id);

  // 开始新batch
  void start_new_batch(int batch_id);

  // 最终化并写入文件
  void finalize(const std::string& filename = "fission_matrix.h5");

  // 伴随源迭代计算（最终化时调用，执行更多迭代）
  // initial_guess: "uniform" 均匀分布, "forward" 正向源分布
  // max_iterations: 最大迭代次数
  // tolerance: 收敛容差
  void compute_adjoint_source(const std::string& initial_guess = "uniform",
    int max_iterations = 1, double tolerance = 1.0e-6);

  // 访问伴随源计算结果
  bool is_adjoint_computed() const { return adjoint_computed_; }
  int get_adjoint_nonzero_cells() const
  {
    return std::count_if(adjoint_source_.begin(), adjoint_source_.end(),
      [](double x) { return x > 1e-10; });
  }

  // 获取标量伴随源数据（用于共轭通量计算，向后兼容）
  const vector<double>& get_adjoint_source() const { return adjoint_source_; }

  // 获取能量分辨伴随源 I*(source_state)；size = n_cells * n_source_groups
  // 若 n_source_groups == 1，与 get_adjoint_source() 相同
  const vector<double>& get_adjoint_source_grouped() const
  {
    return adjoint_source_grouped_;
  }

  // 源能群数（= len(energy_edges)-1，最少为 1）
  int n_source_groups() const { return n_source_groups_; }

  // 源能群边界（与 GreenFunctionMesh 共享）
  const std::vector<double>& source_energy_edges() const
  {
    return source_energy_edges_;
  }

  // source_state 辅助：从 source_state 解码 cell/group
  int state_to_cell(int source_state) const
  {
    return source_state / n_source_groups_;
  }
  int state_to_group(int source_state) const
  {
    return source_state % n_source_groups_;
  }
  int cell_group_to_state(int cell, int group) const
  {
    return cell * n_source_groups_ + group;
  }

  // 网格信息
  const std::array<int, 3>& shape() const { return grid_->shape(); }
  const std::array<double, 3>& origin() const { return grid_->origin(); }
  double pitch() const { return grid_->pitch(); }
  size_t n_cells() const { return grid_->n_cells(); }

private:
  // 将3D位置转换为线性索引
  int position_to_index(const Position& r) const;

  // 由出生能量或 MG 群号确定源能群索引
  int determine_source_group(double energy_eV, int mg_group) const;

  // 稀疏矩阵存储 (COO格式 - Coordinate format)
  // 仅存储非零元素: (row, col, value)
  struct SparseEntry {
    size_t row;   // 源状态索引 source_state = cell * n_source_groups + g
    size_t col;   // 裂变单元索引
    double value; // 裂变中子数
  };

  // 使用map存储当前batch的稀疏数据
  // key = source_state * n_cells + fission_cell
  std::unordered_map<size_t, double> current_batch_sparse_;

  // 累积的稀疏矩阵数据
  std::unordered_map<size_t, double> fission_matrix_sparse_;

  // 记录每个源粒子的出生状态 (source_state = cell * n_source_groups + g)
  std::unordered_map<int64_t, int> source_birth_states_;

  // 每个源状态产生的粒子总数（用于归一化），size = n_source_states
  vector<double> source_counts_;               // 累积
  vector<double> current_batch_source_counts_; // 当前 batch

  // 统一网格配置
  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_; // 保存上边界用于输出
  double inv_pitch_;

  // 源能群配置
  int n_source_groups_ {1};                 // 源中子出生能群数
  std::vector<double> source_energy_edges_; // 源能群边界

  // Batch管理
  int current_batch_id_;
  int max_batches_;
  int n_realizations_; // 已完成的batch数

  // 统计信息
  std::atomic<uint64_t> total_fissions_ {0};
  std::atomic<uint64_t> total_sources_ {0};

  // 伴随源相关
  // adjoint_source_: 标量（按 cell，向后兼容），size = n_cells
  // adjoint_source_grouped_: 能量分辨，size = n_source_states = n_cells *
  // n_source_groups
  vector<double> adjoint_source_;             // I*(cell)，向后兼容
  vector<double> adjoint_source_grouped_;     // I*(source_state)，新增
  vector<double> adjoint_source_batch_;       // 使用每个batch的FM逐步迭代
  vector<double> adjoint_source_accumulated_; // 使用累积的FM逐步迭代
  vector<double> forward_source_;             // 正向源分布 S
  bool adjoint_computed_;                     // 是否已计算伴随源
  int adjoint_iterations_;                    // 实际迭代次数
  double keff_reference_ {1.0};               // 归一化所使用的主计算keff

  // batch级伴随源迭代控制
  bool enable_batch_adjoint_;      // 是否启用每batch迭代
  int adjoint_max_iter_per_batch_; // 每batch迭代次数
  double adjoint_tolerance_;       // 收敛容差
  int adjoint_start_batch_;        // 从哪个batch开始统计FM和计算伴随源

  // 线程安全
  mutable std::mutex data_mutex_;

private:
  // 内部伴随源迭代函数（被batch和finalize调用）
  // iterations: 本次迭代的次数
  // verbose: 是否输出详细信息
  // use_current_batch: true=仅用当前batch的矩阵, false=用累积矩阵
  void perform_adjoint_iteration(
    int iterations, bool verbose = false, bool use_current_batch = false);
};

} // namespace openmc

#endif // OPENMC_FISSION_MATRIX_H
