#ifndef OPENMC_GREENFUNCTION_MESH_H
#define OPENMC_GREENFUNCTION_MESH_H

#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

class GreenFunctionMesh {
public:
  // 构造函数：根据模型边界和分辨率创建网格
  GreenFunctionMesh(double resolution = 1.0); // 默认1cm分辨率
  
  // 根据位置添加贡献
  void accumulate(const Position& r, double contribution);
  
  // 获取网格数据
  const vector<double>& data() const { return data_; }
  
  // 网格尺寸信息
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }

private:
  vector<double> data_;          // 网格数据存储
  std::array<int, 3> shape_;     // 网格尺寸 [nx, ny, nz]
  std::array<double, 3> origin_;  // 网格原点
  double pitch_;                 // 网格间距
};

} // namespace openmc

#endif // OPENMC_GREENFUNCTION_MESH_H