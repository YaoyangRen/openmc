#ifndef OPENMC_GREENFUNCTION_MESH_H
#define OPENMC_GREENFUNCTION_MESH_H

#include "hdf5.h"
#include "openmc/array.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <atomic>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class GreenFunctionMesh {
public:
  // 构造函数，初始化格林函数网格
  explicit GreenFunctionMesh(double resolution, int max_batches);

  // 为特定源粒子累积贡献
  void accumulate(
    const Position& r, double contribution, int64_t source_particle_id);

  // 开始新batch
  void start_new_batch(int batch_id);

  // 写入文件
  void finalize_greenfunction_mesh(const int batch_id);

  // 获取特定源粒子的数据
  const vector<double>& get_particle_data(int64_t source_particle_id) const;

  // 网格尺寸信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }

private:
  // 每个源粒子的格林函数矩阵数据
  // key: source_particle_id, value: 该粒子的格林函数数据
  std::unordered_map<int64_t, vector<double>> particle_green_functions_;

  // 当前batch中每个源粒子的数据
  std::unordered_map<int64_t, vector<double>> current_batch_particle_data_;

  // 累积所有源粒子的数据（总的格林函数）
  vector<double> cumulative_data_;

  std::array<int, 3> shape_;     // 网格的形状（每个维度的单元数）
  std::array<double, 3> origin_; // 网格的原点位置
  double pitch_;                 // 网格单元的边长
  double inv_pitch_;             // 网格单元边长的倒数
  int current_batch_id_;         // 当前处理的batch ID
  int max_batches_;              // 最大batch数量
  size_t spatial_size_;          // 统计信息
  std::atomic<uint64_t> dropped_contributions_ {0};
  std::atomic<uint64_t> total_contributions_ {0};
  
  // 线程同步
  mutable std::mutex data_mutex_;  // 保护current_batch_particle_data_的访问
};

} // namespace openmc

#endif // OPENMC_GREENFUNCTION_MESH_H