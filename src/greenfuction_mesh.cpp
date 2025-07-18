#include "openmc/greenfunction_mesh.h"
#include "openmc/geometry.h"

namespace openmc {

GreenFunctionMesh::GreenFunctionMesh(double resolution) 
  : pitch_(resolution)
{
  // TODO 获取模型边界 这块的参数都只是概念，没有指定
  BoundingBox bbox = model::geometry->bounding_box();
  
  // 设置原点
  origin_ = {bbox.xmin, bbox.ymin, bbox.zmin};
  
  // 计算网格尺寸
  double dx = bbox.xmax - bbox.xmin;
  double dy = bbox.ymax - bbox.ymin;
  double dz = bbox.zmax - bbox.zmin;
  
  shape_[0] = std::ceil(dx / pitch_) + 1;
  shape_[1] = std::ceil(dy / pitch_) + 1;
  shape_[2] = std::ceil(dz / pitch_) + 1;
  
  // 初始化数据数组
  data_.resize(shape_[0] * shape_[1] * shape_[2], 0.0);
}

void GreenFunctionMesh::accumulate(const Position& r, double contribution)
{
  // 计算网格索引
  int i = (r.x - origin_[0]) / pitch_;
  int j = (r.y - origin_[1]) / pitch_;
  int k = (r.z - origin_[2]) / pitch_;
  
  // 边界检查
  if (i < 0 || i >= shape_[0] || 
      j < 0 || j >= shape_[1] ||
      k < 0 || k >= shape_[2]) return;
  
  // 计算一维索引
  int idx = i * shape_[1] * shape_[2] + j * shape_[2] + k;
  
  // 累加贡献
  data_[idx] += contribution;
}

} // namespace openmc