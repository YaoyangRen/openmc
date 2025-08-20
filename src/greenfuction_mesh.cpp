#include "openmc/greenfunction_mesh.h"

#include <algorithm>
#include <cstdint> // for int64_t
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/vector.h"

namespace openmc {

GreenFunctionMesh::GreenFunctionMesh(double resolution)
  : pitch_(resolution), inv_pitch_(1.0 / resolution)
{
  // 调用 C API 获取全局边界
  double llc[3]; // lower left corner
  double urc[3]; // upper right corner
  int err = openmc_global_bounding_box(llc, urc);
  if (err != 0) {
    throw std::runtime_error("Failed to get global bounding box.");
  }

  // 设置原点
  origin_ = {llc[0], llc[1], llc[2]};

  // 计算网格尺寸
  double dx = urc[0] - llc[0];
  double dy = urc[1] - llc[1];
  double dz = urc[2] - llc[2];

  // 确保每个维度至少有一个网格单元
  shape_[0] = std::max(1, static_cast<int>(std::ceil(dx / pitch_)) + 1);
  shape_[1] = std::max(1, static_cast<int>(std::ceil(dy / pitch_)) + 1);
  shape_[2] = std::max(1, static_cast<int>(std::ceil(dz / pitch_)) + 1);

  // 初始化数据数组
  data_.resize(shape_[0] * shape_[1] * shape_[2], 0.0);
}

void GreenFunctionMesh::accumulate(const Position& r, double contribution)
{
  constexpr double eps = 1.0e-10; // 添加容差值

  // 计算网格索引，使用预计算的inv_pitch_提高性能
  int ix = static_cast<int>(std::floor((r.x - origin_[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin_[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin_[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix >= 0 && ix < shape_[0] && iy >= 0 && iy < shape_[1] && iz >= 0 &&
      iz < shape_[2]) {
    // 累积贡献值（线程安全）
    int index = ix + shape_[0] * (iy + shape_[1] * iz);
#pragma omp atomic
    data_[index] += contribution;

    // 记录总贡献数
    total_contributions_++;
  } else {
    // 记录被丢弃的贡献数
    dropped_contributions_++;
  }
}

void GreenFunctionMesh::finalize_greenfunction_mesh(const std::string& filename)
{
  // 创建HDF5文件
  hid_t green_function_file_id = file_open(filename, 'w');

  // 写入文件头部信息
  write_attribute(green_function_file_id, "filetype", "green_function_mesh");
  write_attribute(green_function_file_id, "version", "1.0");
  // write_attribute(green_function_file_id, "date_and_time", time_stamp());

  // 写入网格参数
  write_attribute(green_function_file_id, "pitch", pitch_);
  write_dataset(green_function_file_id, "origin", origin_);
  write_dataset(green_function_file_id, "shape", shape_);

  // 写入数据
  write_dataset(green_function_file_id, "green_function_data", data_);

  // 关闭文件
  file_close(green_function_file_id);

  // 清除data_以释放内存
  data_.clear();
  data_.shrink_to_fit();
}

} // namespace openmc