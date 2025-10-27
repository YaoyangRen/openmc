#ifndef OPENMC_FISSION_MATRIX_H
#define OPENMC_FISSION_MATRIX_H

#include "hdf5.h"
#include "openmc/array.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <atomic>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class FissionMatrix {
public:
  // 构造函数
  explicit FissionMatrix(double resolution, int max_batches,
    bool auto_bounds = true,
    const std::array<double, 3>& manual_lower = {0.0, 0.0, 0.0},
    const std::array<double, 3>& manual_upper = {10.0, 10.0, 10.0});

  // 记录源中子的产生位置
  void record_source_birth(const Position& r, int64_t source_particle_id);

  // 记录裂变事件（裂变位置和产生的中子数）
  void record_fission_event(
    const Position& r, double nu_fission, int64_t source_particle_id);

  // 开始新batch
  void start_new_batch(int batch_id);

  // 最终化并写入文件
  void finalize(const std::string& filename = "fission_matrix.h5");

  // 网格信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }
  size_t n_cells() const { return n_cells_; }

private:
  // 将3D位置转换为线性索引
  int position_to_index(const Position& r) const;

  // 裂变矩阵数据: F[i][j] = 源在i产生, 在j裂变的中子数
  vector<vector<double>> fission_matrix_;

  // 当前batch的裂变矩阵
  vector<vector<double>> current_batch_matrix_;

  // 记录每个源粒子的出生位置
  std::unordered_map<int64_t, int> source_birth_cells_;

  // 每个源区域产生的粒子总数（用于归一化）
  vector<double> source_counts_;
  vector<double> current_batch_source_counts_;

  // 网格参数
  std::array<int, 3> shape_;
  std::array<double, 3> origin_;
  std::array<double, 3> upper_bound_; // 保存上边界用于输出
  double pitch_;
  double inv_pitch_;
  size_t n_cells_; // 总单元数

  // Batch管理
  int current_batch_id_;
  int max_batches_;
  int n_realizations_; // 已完成的batch数

  // 统计信息
  std::atomic<uint64_t> total_fissions_ {0};
  std::atomic<uint64_t> total_sources_ {0};

  // 线程安全
  mutable std::mutex data_mutex_;
};

} // namespace openmc

#endif // OPENMC_FISSION_MATRIX_H
