#ifndef OPENMC_MESH_INIT_H
#define OPENMC_MESH_INIT_H

#include "openmc/array.h"
#include "openmc/position.h"
#include <memory>
#include <string>

namespace openmc {

//! \class SharedMeshGrid
//! \brief 裂变矩阵和传递函数共享的统一网格
//!
//! 该类提供统一的空间网格划分，确保裂变矩阵和传递函数使用完全相同的网格参数
//! 避免网格不一致导致的伴随源计算错误
//!
class SharedMeshGrid {
public:
  // 静态工厂方法：创建统一网格
  // resolution: 网格分辨率 (cm)
  // auto_bounds: 是否自动从几何体获取边界
  // manual_lower: 手动指定的下边界 [x, y, z]
  // manual_upper: 手动指定的上边界 [x, y, z]
  static std::shared_ptr<SharedMeshGrid> create(double resolution,
    bool auto_bounds = true,
    const std::array<double, 3>& manual_lower = {0.0, 0.0, 0.0},
    const std::array<double, 3>& manual_upper = {10.0, 10.0, 10.0});

  // Getter 方法
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  const std::array<double, 3>& upper_bound() const { return upper_bound_; }
  double pitch() const { return pitch_; }
  double inv_pitch() const { return inv_pitch_; }
  size_t n_cells() const { return n_cells_; }

  // 验证网格配置
  void validate() const;

  // 输出网格信息
  void print_info(const std::string& label = "Shared Mesh Grid") const;

private:
  // 私有构造函数（通过 create() 工厂方法创建）
  SharedMeshGrid(double resolution, bool auto_bounds,
    const std::array<double, 3>& manual_lower,
    const std::array<double, 3>& manual_upper);

  // 初始化网格（从根宇宙或手动边界计算网格参数）
  void initialize();

  // 网格参数
  std::array<int, 3> shape_ {0, 0, 0};                // 网格形状 [nx, ny, nz]
  std::array<double, 3> origin_ {0.0, 0.0, 0.0};      // 原点坐标
  std::array<double, 3> upper_bound_ {0.0, 0.0, 0.0}; // 上边界坐标
  double pitch_ {0.0};                                // 网格分辨率 (cm)
  double inv_pitch_ {0.0};                            // 分辨率的倒数
  size_t n_cells_ {0};                                // 总单元数

  // 配置参数
  bool auto_bounds_ {true};            // 是否自动获取边界
  std::array<double, 3> manual_lower_; // 手动下边界
  std::array<double, 3> manual_upper_; // 手动上边界
};

} // namespace openmc

#endif // OPENMC_MESH_INIT_H
