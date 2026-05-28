#include "openmc/mesh_init.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/settings.h"
#include "openmc/universe.h"

namespace openmc {

// 静态工厂方法实现
std::shared_ptr<SharedMeshGrid> SharedMeshGrid::create(double resolution,
  bool auto_bounds, const std::array<double, 3>& manual_lower,
  const std::array<double, 3>& manual_upper, bool has_manual_bounds)
{
  if (resolution <= 0.0) {
    fatal_error("SharedMeshGrid resolution must be positive.");
  }
  auto grid = std::shared_ptr<SharedMeshGrid>(
    new SharedMeshGrid(
      resolution, auto_bounds, manual_lower, manual_upper, has_manual_bounds));
  grid->initialize();
  return grid;
}

std::shared_ptr<SharedMeshGrid> SharedMeshGrid::create_from_settings()
{
  return create(settings::kinetics_mesh_pitch, settings::kinetics_mesh_auto_bounds,
    settings::kinetics_mesh_lower_left, settings::kinetics_mesh_upper_right,
    settings::kinetics_mesh_has_bounds);
}

// 私有构造函数
SharedMeshGrid::SharedMeshGrid(double resolution, bool auto_bounds,
  const std::array<double, 3>& manual_lower,
  const std::array<double, 3>& manual_upper, bool has_manual_bounds)
  : has_manual_bounds_(has_manual_bounds), pitch_(resolution),
    inv_pitch_(1.0 / resolution), auto_bounds_(auto_bounds),
    manual_lower_(manual_lower), manual_upper_(manual_upper)
{}

// 初始化网格
void SharedMeshGrid::initialize()
{
  double llc[3];
  double urc[3];

  if (auto_bounds_) {
    // 自动从根宇宙获取边界
    auto bbox = model::universes.at(model::root_universe)->bounding_box();

    llc[0] = bbox.xmin;
    llc[1] = bbox.ymin;
    llc[2] = bbox.zmin;

    urc[0] = bbox.xmax;
    urc[1] = bbox.ymax;
    urc[2] = bbox.zmax;

    // 检查无限边界
    bool has_infinite = false;
    for (int i = 0; i < 3; ++i) {
      if (llc[i] <= -INFTY || urc[i] >= INFTY) {
        has_infinite = true;
        break;
      }
    }

    if (has_infinite) {
      if (!has_manual_bounds_) {
        fatal_error("Kinetics mesh automatic bounds found an infinite geometry "
                    "extent. Specify <kinetics_mesh>/<lower_left> and "
                    "<upper_right> in settings.xml, or set "
                    "<auto_bounds>false.");
      }
      // 使用手动边界
      bounds_source_ =
        "manual bounds from settings.xml (geometry XML bounds infinite)";
      std::copy(manual_lower_.begin(), manual_lower_.end(), llc);
      std::copy(manual_upper_.begin(), manual_upper_.end(), urc);
    } else {
      bounds_source_ =
        "XML geometry bounding box (read from geometry/model XML)";
      // 添加5%边距
      double margin = pitch_ * 0.05;
      for (int i = 0; i < 3; ++i) {
        llc[i] -= margin;
        urc[i] += margin;
      }
    }
  } else {
    if (!has_manual_bounds_) {
      fatal_error("Kinetics mesh manual bounds were requested, but "
                  "<kinetics_mesh>/<lower_left> and <upper_right> were not "
                  "specified.");
    }
    // 使用手动边界
    bounds_source_ = "manual bounds from settings.xml (auto_bounds=false)";
    std::copy(manual_lower_.begin(), manual_lower_.end(), llc);
    std::copy(manual_upper_.begin(), manual_upper_.end(), urc);
  }

  origin_ = {llc[0], llc[1], llc[2]};
  upper_bound_ = {urc[0], urc[1], urc[2]};

  // 计算网格尺寸
  double dx = urc[0] - llc[0];
  double dy = urc[1] - llc[1];
  double dz = urc[2] - llc[2];

  shape_[0] = std::max(1, static_cast<int>(std::ceil(dx / pitch_)) + 1);
  shape_[1] = std::max(1, static_cast<int>(std::ceil(dy / pitch_)) + 1);
  shape_[2] = std::max(1, static_cast<int>(std::ceil(dz / pitch_)) + 1);

  n_cells_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  print_info();
}

// 验证网格配置
void SharedMeshGrid::validate() const
{
  if (n_cells_ == 0) {
    std::cerr << "[SharedMeshGrid] 错误：网格未初始化\n";
    return;
  }

  if (pitch_ <= 0) {
    std::cerr << "[SharedMeshGrid] 错误：无效的网格分辨率 " << pitch_
              << " cm\n";
    return;
  }

  for (int i = 0; i < 3; ++i) {
    if (shape_[i] <= 0) {
      std::cerr << "[SharedMeshGrid] 错误：无效的网格形状 [" << shape_[0]
                << ", " << shape_[1] << ", " << shape_[2] << "]\n";
      return;
    }
  }
}

// 输出网格信息
void SharedMeshGrid::print_info(const std::string& label) const
{
  std::cout << "\n[" << label << "] 几何网格划分" << std::endl;
  std::cout << "  shape: " << shape_[0] << " x " << shape_[1] << " x "
            << shape_[2] << ", pitch: " << pitch_ << " cm, cells: "
            << n_cells_ << std::endl;
  std::cout << "  bounds_source: " << bounds_source_ << std::endl;
  std::cout << "  lower_left: [" << origin_[0] << ", " << origin_[1] << ", "
            << origin_[2] << "] cm, upper_right: [" << upper_bound_[0]
            << ", " << upper_bound_[1] << ", " << upper_bound_[2]
            << "] cm" << std::endl;
}

} // namespace openmc
