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
  explicit GreenFunctionMesh(double resolution = 1.0);

  // 根据位置添加贡献
  void accumulate(const Position& r, double contribution);

  // 写入文件
  void finalize_greenfunction_mesh(const std::string& filename);

  // 获取网格数据
  const vector<double>& data() const { return data_; }

  // 网格尺寸信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }

  // 获取统计信息
  uint64_t get_dropped_contributions() const { return dropped_contributions_; }
  uint64_t get_total_contributions() const { return total_contributions_; }

private:
  vector<double> data_;          // 网格数据存储
  std::array<int, 3> shape_;     // 网格尺寸 [nx, ny, nz]
  std::array<double, 3> origin_; // 网格原点
  double pitch_;                 // 网格间距
  double inv_pitch_;             // 网格间距的倒数，用于优化计算
  std::atomic<uint64_t> dropped_contributions_ {0}; // 记录被丢弃的贡献数
  std::atomic<uint64_t> total_contributions_ {0};   // 记录总的贡献数
};

} // namespace openmc

#endif // OPENMC_GREENFUNCTION_MESH_H