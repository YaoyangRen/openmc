#include "openmc/flux_mesh.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh_init.h"
#include "openmc/message_passing.h"

namespace openmc {

//==============================================================================
// FluxMesh 实现
//==============================================================================

FluxMesh::FluxMesh(std::shared_ptr<SharedMeshGrid> grid)
  : grid_(grid), n_accumulated_batches_(0)
{
  if (!grid_) {
    fatal_error("FluxMesh: 共享网格指针为空!");
  }

  // 为每个 OpenMP 线程分配独立的缓冲区
#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif
  thread_flux_.resize(n_threads);

  std::cout << "通量网格已初始化:" << std::endl;
  std::cout << "  网格尺寸: " << grid_->shape()[0] << " x " << grid_->shape()[1]
            << " x " << grid_->shape()[2] << std::endl;
  std::cout << "  总单元数: " << grid_->n_cells() << std::endl;
  std::cout << "  边界范围: [" << grid_->origin()[0] << ", "
            << grid_->origin()[1] << ", " << grid_->origin()[2] << "] 到 ["
            << grid_->upper_bound()[0] << ", " << grid_->upper_bound()[1]
            << ", " << grid_->upper_bound()[2] << "]" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

//------------------------------------------------------------------------------

void FluxMesh::accumulate(
  const std::array<double, 3>& position, double weight, double distance)
{
  // 获取网格单元索引
  int cell_index = position_to_index(position);

  // 如果位置在网格范围内
  if (cell_index >= 0) {
    // 径迹长度估计器: Φ = w × d
    double contribution = weight * distance;

    // 获取当前线程ID并累积到线程局部缓冲区(无锁)
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    thread_flux_[thread_id][cell_index] += contribution;
  }
}

//------------------------------------------------------------------------------

void FluxMesh::end_batch(int batch)
{
  // 首先合并所有线程的局部数据
  merge_thread_local_data();

  // 跨进程归约当前批次的通量数据
#ifdef OPENMC_MPI
  // 收集所有进程的稀疏数据
  // 注意: 这里简化处理,实际应该使用MPI归约操作
  // 对于稀疏数据,可以先收集所有非零索引,再归约对应的值
#endif

  // 将当前批次的结果累积到统计量中
  for (const auto& [cell_index, flux_value] : batch_flux_) {
    flux_sum_[cell_index] += flux_value;
    flux_sum_sq_[cell_index] += flux_value * flux_value;
  }

  // 清空当前批次数据
  batch_flux_.clear();

  // 增加已累积的批次计数
  n_accumulated_batches_++;

  // 输出进度信息
  // if (batch % 10 == 0) {
  //   std::cout << "通量网格: 批次 " << batch << " 已处理 ("
  //             << n_accumulated_batches_ << " 批次已累积)" << std::endl;
  // }
}

//------------------------------------------------------------------------------

void FluxMesh::merge_thread_local_data()
{
  // 合并所有线程的局部数据到 batch_flux_
  for (auto& thread_map : thread_flux_) {
    for (const auto& [cell_index, flux_value] : thread_map) {
      batch_flux_[cell_index] += flux_value;
    }
    // 清空线程局部数据,准备下一批次
    thread_map.clear();
  }
}

//------------------------------------------------------------------------------

void FluxMesh::finalize(int n_batches)
{
  std::cout << "通量网格最终化: 共处理 " << n_accumulated_batches_ << " 批次"
            << std::endl;
  std::cout << "  非零通量网格单元数: " << flux_sum_.size() << std::endl;

  // 输出到 HDF5 文件
  write_hdf5("flux_mesh.h5", n_batches);

  std::cout << "通量分布已保存到 flux_mesh.h5" << std::endl;
}

//------------------------------------------------------------------------------

void FluxMesh::reset()
{
  batch_flux_.clear();
  for (auto& thread_map : thread_flux_) {
    thread_map.clear();
  }
  flux_sum_.clear();
  flux_sum_sq_.clear();
  n_accumulated_batches_ = 0;
}

//------------------------------------------------------------------------------

double FluxMesh::get_flux(int cell_index) const
{
  auto it = flux_sum_.find(cell_index);
  if (it != flux_sum_.end() && n_accumulated_batches_ > 0) {
    return it->second / n_accumulated_batches_;
  }
  return 0.0;
}

//------------------------------------------------------------------------------

int FluxMesh::position_to_index(const std::array<double, 3>& position) const
{
  const auto& origin = grid_->origin();
  const auto& upper_bound = grid_->upper_bound();
  const auto& shape = grid_->shape();
  double inv_pitch = grid_->inv_pitch();

  // 检查位置是否在网格范围内
  for (int i = 0; i < 3; ++i) {
    if (position[i] < origin[i] || position[i] > upper_bound[i]) {
      return -1; // 超出范围
    }
  }

  // 计算三维网格索引
  int ix = static_cast<int>((position[0] - origin[0]) * inv_pitch);
  int iy = static_cast<int>((position[1] - origin[1]) * inv_pitch);
  int iz = static_cast<int>((position[2] - origin[2]) * inv_pitch);

  // 边界检查(处理舍入误差)
  ix = std::min(ix, shape[0] - 1);
  iy = std::min(iy, shape[1] - 1);
  iz = std::min(iz, shape[2] - 1);

  // 转换为一维索引
  return grid_to_index(ix, iy, iz);
}

//------------------------------------------------------------------------------

std::array<int, 3> FluxMesh::index_to_grid(int index) const
{
  const auto& shape = grid_->shape();
  std::array<int, 3> grid_idx;
  grid_idx[2] = index / (shape[0] * shape[1]); // iz
  int remainder = index % (shape[0] * shape[1]);
  grid_idx[1] = remainder / shape[0]; // iy
  grid_idx[0] = remainder % shape[0]; // ix
  return grid_idx;
}

//------------------------------------------------------------------------------

int FluxMesh::grid_to_index(int ix, int iy, int iz) const
{
  const auto& shape = grid_->shape();
  return iz * shape[0] * shape[1] + iy * shape[0] + ix;
}

//------------------------------------------------------------------------------

void FluxMesh::write_hdf5(const std::string& filename, int n_batches) const
{
  // 只有主进程写入文件
#ifdef OPENMC_MPI
  if (mpi::rank != 0)
    return;
#endif

  // 创建 HDF5 文件
  hid_t file_id = file_open(filename, 'w');

  // 获取网格参数
  const auto& shape = grid_->shape();
  const auto& origin = grid_->origin();
  const auto& upper_bound = grid_->upper_bound();
  double pitch = grid_->pitch();
  size_t n_cells = grid_->n_cells();

  // 写入网格参数
  std::array<int, 3> grid_shape = {shape[0], shape[1], shape[2]};
  std::array<double, 3> grid_lower = {origin[0], origin[1], origin[2]};
  std::array<double, 3> grid_upper = {
    upper_bound[0], upper_bound[1], upper_bound[2]};
  std::array<double, 3> grid_pitch = {pitch, pitch, pitch};

  write_dataset(file_id, "grid_shape", grid_shape);
  write_dataset(file_id, "grid_lower_left", grid_lower);
  write_dataset(file_id, "grid_upper_right", grid_upper);
  write_dataset(file_id, "grid_pitch", grid_pitch);

  int total_cells = static_cast<int>(n_cells);
  write_dataset(file_id, "n_cells", total_cells);
  write_dataset(file_id, "n_batches", n_batches);

  // 准备稀疏数据数组
  std::vector<int> cell_indices;
  std::vector<double> flux_mean;
  std::vector<double> flux_std;

  cell_indices.reserve(flux_sum_.size());
  flux_mean.reserve(flux_sum_.size());
  flux_std.reserve(flux_sum_.size());

  // 计算每个非零单元的均值和标准差
  for (const auto& [cell_index, sum] : flux_sum_) {
    double mean = sum / n_accumulated_batches_;
    double sum_sq = flux_sum_sq_.at(cell_index);
    double variance = 0.0;

    // 计算标准差: σ = sqrt(E[X²] - E[X]²)
    if (n_accumulated_batches_ > 1) {
      double mean_sq = sum_sq / n_accumulated_batches_;
      variance = mean_sq - mean * mean;
      // 防止数值误差导致负方差
      variance = std::max(variance, 0.0);
    }
    double std_dev = std::sqrt(variance);

    cell_indices.push_back(cell_index);
    flux_mean.push_back(mean);
    flux_std.push_back(std_dev);
  }

  // 写入稀疏数据
  write_dataset(file_id, "cell_indices", cell_indices);
  write_dataset(file_id, "flux_mean", flux_mean);
  write_dataset(file_id, "flux_std", flux_std);

  // 创建密集格式的通量数组(用于可视化)
  std::vector<double> flux_dense(n_cells, 0.0);
  std::vector<double> flux_error_dense(n_cells, 0.0);

  for (size_t i = 0; i < cell_indices.size(); ++i) {
    int idx = cell_indices[i];
    if (idx >= 0 && idx < static_cast<int>(n_cells)) {
      flux_dense[idx] = flux_mean[i];
      flux_error_dense[idx] = flux_std[i];
    }
  }

  // 写入密集格式数据
  write_dataset(file_id, "flux_mean_dense", flux_dense);
  write_dataset(file_id, "flux_std_dense", flux_error_dense);

  // 关闭文件
  file_close(file_id);
}

} // namespace openmc
