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
  explicit FissionMatrix(std::shared_ptr<SharedMeshGrid> grid, int max_batches);

  // 记录源中子的产生位置
  void record_source_birth(const Position& r, int64_t source_particle_id);

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

  // 启用/禁用每个batch后的伴随源迭代
  // enable: 是否启用batch级伴随源迭代
  // iterations_per_batch: 每个batch执行的迭代次数
  // tolerance: 收敛容差
  // start_batch: 从哪个batch开始统计（默认从第5个batch开始）
  void enable_batch_adjoint_iteration(bool enable = true,
    int iterations_per_batch = 10, double tolerance = 1.0e-6,
    int start_batch = 5);

  // 计算伴随源分布（用于最后一个非活跃batch）
  void compute_batch_adjoint();

  // 访问伴随源计算结果
  bool is_adjoint_computed() const { return adjoint_computed_; }
  int get_adjoint_nonzero_cells() const
  {
    return std::count_if(adjoint_source_.begin(), adjoint_source_.end(),
      [](double x) { return x > 1e-10; });
  }

  // 获取伴随源数据（用于共轭通量计算）
  const vector<double>& get_adjoint_source() const { return adjoint_source_; }

  // 网格信息
  const std::array<int, 3>& shape() const { return grid_->shape(); }
  const std::array<double, 3>& origin() const { return grid_->origin(); }
  double pitch() const { return grid_->pitch(); }
  size_t n_cells() const { return grid_->n_cells(); }

private:
  // 将3D位置转换为线性索引
  int position_to_index(const Position& r) const;

  // 稀疏矩阵存储 (COO格式 - Coordinate format)
  // 仅存储非零元素: (row, col, value)
  struct SparseEntry {
    size_t row;   // 源单元索引
    size_t col;   // 裂变单元索引
    double value; // 裂变中子数
  };

  // 使用map存储当前batch的稀疏数据，key = row * n_cells + col
  std::unordered_map<size_t, double> current_batch_sparse_;

  // 累积的稀疏矩阵数据
  std::unordered_map<size_t, double> fission_matrix_sparse_;

  // 记录每个源粒子的出生位置
  std::unordered_map<int64_t, int> source_birth_cells_;

  // 每个源区域产生的粒子总数（用于归一化）
  vector<double> source_counts_;
  vector<double> current_batch_source_counts_;

  // 统一网格配置
  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_; // 保存上边界用于输出
  double inv_pitch_;

  // Batch管理
  int current_batch_id_;
  int max_batches_;
  int n_realizations_; // 已完成的batch数

  // 统计信息
  std::atomic<uint64_t> total_fissions_ {0};
  std::atomic<uint64_t> total_sources_ {0};

  // 伴随源相关
  vector<double> adjoint_source_;       // 伴随源分布 I* (默认，使用累积FM)
  vector<double> adjoint_source_batch_; // 使用每个batch的FM逐步迭代
  vector<double> adjoint_source_accumulated_; // 使用累积的FM逐步迭代
  vector<double> forward_source_;             // 正向源分布 S (用于初始化)
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
