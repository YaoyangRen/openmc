#!/usr/bin/env python3
"""修改 greenfunction_mesh.cpp 使用 SharedMeshGrid"""

with open('d:/OpenMC/openmc/src/greenfunction_mesh.cpp', 'r', encoding='utf-8') as f:
    content = f.read()

# 修改构造函数签名
content = content.replace(
    'GreenFunctionMesh::GreenFunctionMesh(double resolution, int max_batches,\n  bool auto_bounds, const std::array<double, 3>& manual_lower,\n  const std::array<double, 3>& manual_upper)\n  : pitch_(resolution), inv_pitch_(1.0 / resolution), current_batch_id_(-1),\n    max_batches_(max_batches)',
    'GreenFunctionMesh::GreenFunctionMesh(std::shared_ptr<SharedMeshGrid> grid, int max_batches)\n  : grid_(grid), inv_pitch_(1.0 / grid->pitch()), current_batch_id_(-1),\n    max_batches_(max_batches)'
)

# 删除旧的网格初始化代码
import re
old_init = re.compile(r'{\n\n  double llc\[3\];.*?spatial_size_ = static_cast<size_t>\(shape_\[0\]\) \* shape_\[1\] \* shape_\[2\];', re.DOTALL)
new_init = '''{
  // 计算上边界（用于输出）
  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();
  double pitch = grid_->pitch();
  upper_bound_ = {origin[0] + shape[0] * pitch,
                  origin[1] + shape[1] * pitch,
                  origin[2] + shape[2] * pitch};

  spatial_size_ = grid_->n_cells();'''
content = old_init.sub(new_init, content)

# 修改 position_to_cell_index
old_pos_func = '''int GreenFunctionMesh::position_to_cell_index(const Position& r) const
{
  constexpr double eps = 1.0e-10;

  int ix = static_cast<int>(std::floor((r.x - origin_[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin_[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin_[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix < 0 || ix >= shape_[0] || iy < 0 || iy >= shape_[1] || iz < 0 ||
      iz >= shape_[2]) {
    return -1;
  }

  // 计算线性索引：ix + iy*nx + iz*nx*ny
  int index = ix + shape_[0] * (iy + shape_[1] * iz);
  return index;
}'''

new_pos_func = '''int GreenFunctionMesh::position_to_cell_index(const Position& r) const
{
  constexpr double eps = 1.0e-10;

  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();

  int ix = static_cast<int>(std::floor((r.x - origin[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix < 0 || ix >= shape[0] || iy < 0 || iy >= shape[1] || iz < 0 ||
      iz >= shape[2]) {
    return -1;
  }

  // 计算线性索引：ix + iy*nx + iz*nx*ny
  return ix + iy * shape[0] + iz * shape[0] * shape[1];
}'''
content = content.replace(old_pos_func, new_pos_func)

# 替换其他 shape_, pitch_ 引用（不替换 spatial_size_）
content = content.replace('shape_[', 'grid_->shape()[')
content = content.replace('pitch_', 'grid_->pitch()')

# 保存
with open('d:/OpenMC/openmc/src/greenfunction_mesh.cpp', 'w', encoding='utf-8') as f:
    f.write(content)

print("修改完成！")
