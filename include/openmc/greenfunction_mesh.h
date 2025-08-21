#ifndef OPENMC_GREENFUNCTION_MESH_H
#define OPENMC_GREENFUNCTION_MESH_H

#include "hdf5.h"

#include "openmc/array.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <atomic>
#include <string>

namespace openmc {

class GreenFunctionMesh {
public:
  // 构造函数：根据模型边界和分辨率创建网格
  explicit GreenFunctionMesh(double resolution, int max_batches);

  // 根据位置添加贡献
  void accumulate(const Position& r, double contribution);

  // 开始新batch
  void start_new_batch(int batch_id);

  // 写入文件
  void finalize_greenfunction_mesh(const int batch_id);

  // 获取网格数据
  const vector<double>& data() const { return current_batch_data_; }

  // 网格尺寸信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }

private:
  vector<vector<double>>
    batch_data_; // 每个batch的网格数据 [batch][spatial_index]
  vector<double> current_batch_data_; // 当前batch的数据缓存
  std::array<int, 3> shape_;          // 网格尺寸 [nx, ny, nz]
  std::array<double, 3> origin_;      // 网格原点
  double pitch_;                      // 网格间距
  double inv_pitch_;                  // 网格间距的倒数
  int current_batch_id_;              // 当前batch ID
  int max_batches_;                   // 最大batch数量
  std::atomic<uint64_t> dropped_contributions_ {0};
  std::atomic<uint64_t> total_contributions_ {0};
};

} // namespace openmc

#endif // OPENMC_GREENFUNCTION_MESH_H