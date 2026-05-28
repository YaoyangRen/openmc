#include "openmc/flux_mesh.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
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

FluxMesh::FluxMesh(
  std::shared_ptr<SharedMeshGrid> grid, std::vector<double> energy_edges)
  : grid_(grid), n_accumulated_batches_(0),
    energy_edges_(std::move(energy_edges))
{
  if (!grid_) {
    fatal_error("FluxMesh: 共享网格指针为空!");
  }

  if (energy_edges_.size() >= 2) {
    if (!std::is_sorted(energy_edges_.begin(), energy_edges_.end())) {
      fatal_error(
        "FluxMesh: energy group edges must be sorted in ascending order.");
    }
    n_groups_ = static_cast<int>(energy_edges_.size()) - 1;
  } else {
    energy_edges_.clear();
    n_groups_ = 1;
  }

  // 为每个 OpenMP 线程分配独立的缓冲区
#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif
  thread_flux_group_.resize(n_threads);
}

//------------------------------------------------------------------------------

void FluxMesh::accumulate(const std::array<double, 3>& position, double weight,
  double distance, double energy_eV, int mg_group)
{
  // 获取网格单元索引
  int cell_index = position_to_index(position);

  // 如果位置在网格范围内
  if (cell_index >= 0) {
    int group = determine_group(energy_eV, mg_group);
    if (group < 0) {
      fatal_error("FluxMesh: 无法确定通量统计的能群,请提供有效能量或群索引.");
    }

    // 径迹长度估计器: Φ = w × d
    double contribution = weight * distance;

    // 获取当前线程ID并累积到线程局部缓冲区(无锁)
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    auto& cell_vector = thread_flux_group_[thread_id][cell_index];
    if (cell_vector.empty()) {
      cell_vector = make_zero_group_vector();
    }
    cell_vector[group] += contribution;
  }
}

//------------------------------------------------------------------------------

void FluxMesh::accumulate_track(const std::array<double, 3>& start_position,
  const std::array<double, 3>& direction, double weight, double distance,
  double energy_eV, int mg_group)
{
  if (distance <= 0.0 || weight == 0.0) {
    return;
  }

  int group = determine_group(energy_eV, mg_group);
  if (group < 0) {
    fatal_error("FluxMesh: unable to determine energy group for track-length "
                "flux scoring.");
  }

  const auto& origin = grid_->origin();
  const auto& upper_bound = grid_->upper_bound();
  double pitch = grid_->pitch();
  constexpr double eps = 1.0e-12;

  double t_min = 0.0;
  double t_max = distance;
  for (int axis = 0; axis < 3; ++axis) {
    double x0 = start_position[axis];
    double u = direction[axis];
    if (std::abs(u) < eps) {
      if (x0 < origin[axis] || x0 >= upper_bound[axis]) {
        return;
      }
      continue;
    }

    double t1 = (origin[axis] - x0) / u;
    double t2 = (upper_bound[axis] - x0) / u;
    if (t1 > t2) {
      std::swap(t1, t2);
    }
    t_min = std::max(t_min, t1);
    t_max = std::min(t_max, t2);
    if (t_min >= t_max) {
      return;
    }
  }

  double t = std::max(0.0, t_min);
  while (t < t_max - eps) {
    double probe_t = std::min(t + eps * std::max(1.0, distance), t_max);
    std::array<double, 3> probe {};
    for (int axis = 0; axis < 3; ++axis) {
      probe[axis] = start_position[axis] + probe_t * direction[axis];
    }

    int cell_index = position_to_index(probe);
    if (cell_index < 0) {
      t = probe_t;
      continue;
    }

    auto grid_idx = index_to_grid(cell_index);
    double next_t = t_max;
    for (int axis = 0; axis < 3; ++axis) {
      double u = direction[axis];
      if (std::abs(u) < eps) {
        continue;
      }

      double boundary_coord = 0.0;
      if (u > 0.0) {
        boundary_coord = origin[axis] + (grid_idx[axis] + 1) * pitch;
        boundary_coord = std::min(boundary_coord, upper_bound[axis]);
      } else {
        boundary_coord = origin[axis] + grid_idx[axis] * pitch;
        boundary_coord = std::max(boundary_coord, origin[axis]);
      }

      double boundary_t = (boundary_coord - start_position[axis]) / u;
      if (boundary_t > t + eps) {
        next_t = std::min(next_t, boundary_t);
      }
    }

    if (next_t <= t + eps) {
      next_t = std::min(t_max, t + eps * std::max(1.0, distance));
    }

    double segment = next_t - t;
    if (segment > 0.0) {
      add_contribution(cell_index, group, weight * segment);
    }
    t = next_t;
  }
}

//------------------------------------------------------------------------------

void FluxMesh::add_contribution(
  int cell_index, int group, double contribution)
{
  if (cell_index < 0 || contribution == 0.0) {
    return;
  }

#ifdef _OPENMP
  int thread_id = omp_get_thread_num();
#else
  int thread_id = 0;
#endif
  auto& cell_vector = thread_flux_group_[thread_id][cell_index];
  if (cell_vector.empty()) {
    cell_vector = make_zero_group_vector();
  }
  cell_vector[group] += contribution;
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
  for (const auto& [cell_index, group_values] : batch_flux_group_) {
    double total_value =
      std::accumulate(group_values.begin(), group_values.end(), 0.0);
    flux_sum_[cell_index] += total_value;
    flux_sum_sq_[cell_index] += total_value * total_value;

    auto& group_sum = flux_group_sum_[cell_index];
    if (group_sum.empty()) {
      group_sum = make_zero_group_vector();
    }
    auto& group_sum_sq = flux_group_sum_sq_[cell_index];
    if (group_sum_sq.empty()) {
      group_sum_sq = make_zero_group_vector();
    }

    for (int g = 0; g < n_groups_; ++g) {
      double value = group_values[g];
      group_sum[g] += value;
      group_sum_sq[g] += value * value;
    }
  }

  // 清空当前批次数据
  batch_flux_group_.clear();

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
  // 合并所有线程的局部数据到 batch_flux_group_
  for (auto& thread_map : thread_flux_group_) {
    for (const auto& [cell_index, group_values] : thread_map) {
      auto& batch_vector = batch_flux_group_[cell_index];
      if (batch_vector.empty()) {
        batch_vector = make_zero_group_vector();
      }
      for (int g = 0; g < n_groups_; ++g) {
        batch_vector[g] += group_values[g];
      }
    }
    // 清空线程局部数据,准备下一批次
    thread_map.clear();
  }
}

//------------------------------------------------------------------------------

void FluxMesh::finalize(int n_batches)
{
  std::cout << std::string(70, '=') << std::endl;
  std::cout << "通量网格最终化: 共处理 " << n_accumulated_batches_ << " 批次"
            << std::endl;
  std::cout << "  非零通量网格单元数: " << flux_sum_.size() << std::endl;

  if (n_accumulated_batches_ > 0 && !flux_sum_.empty()) {
    double cell_volume = grid_->pitch() * grid_->pitch() * grid_->pitch();
    double total_flux_mean = 0.0;
    std::vector<double> group_totals(static_cast<size_t>(n_groups_), 0.0);
    std::vector<double> zero_vector = make_zero_group_vector();

    for (const auto& [cell_index, sum] : flux_sum_) {
      double cell_mean = sum / n_accumulated_batches_ / cell_volume;
      total_flux_mean += cell_mean;

      if (n_groups_ > 0) {
        const auto& group_sum = flux_group_sum_.count(cell_index)
                                  ? flux_group_sum_.at(cell_index)
                                  : zero_vector;
        for (int g = 0; g < n_groups_; ++g) {
          group_totals[g] +=
            group_sum[g] / n_accumulated_batches_ / cell_volume;
        }
      }
    }

    std::cout << "  总通量均值(所有单元求和): " << total_flux_mean << std::endl;

    if (n_groups_ > 1) {
      std::cout << "  分群统计:" << std::endl;
      for (int g = 0; g < n_groups_; ++g) {
        double fraction =
          total_flux_mean > 0.0 ? group_totals[g] / total_flux_mean : 0.0;
        if (!energy_edges_.empty() &&
            static_cast<int>(energy_edges_.size()) == n_groups_ + 1) {
          std::cout << "    G" << std::setw(2) << g << " [" << std::scientific
                    << std::setprecision(3) << energy_edges_[g] << ", "
                    << energy_edges_[g + 1] << ") eV: " << group_totals[g]
                    << " (" << std::fixed << std::setprecision(2)
                    << fraction * 100.0 << "%)" << std::endl;
        } else {
          std::cout << "    G" << std::setw(2) << g << ": " << std::scientific
                    << std::setprecision(6) << group_totals[g] << " ("
                    << std::fixed << std::setprecision(2) << fraction * 100.0
                    << "%)" << std::endl;
        }
      }
      std::cout << std::defaultfloat << std::setprecision(6);
    }
  }

  // 输出到 HDF5 文件
  write_hdf5("flux_mesh.h5", n_batches);

  std::cout << "通量分布已保存到 flux_mesh.h5" << std::endl;
}

//------------------------------------------------------------------------------

void FluxMesh::reset()
{
  batch_flux_group_.clear();
  for (auto& thread_map : thread_flux_group_) {
    thread_map.clear();
  }
  flux_group_sum_.clear();
  flux_group_sum_sq_.clear();
  flux_sum_.clear();
  flux_sum_sq_.clear();
  n_accumulated_batches_ = 0;
}

//------------------------------------------------------------------------------

double FluxMesh::get_flux(int cell_index) const
{
  auto it = flux_sum_.find(cell_index);
  if (it != flux_sum_.end() && n_accumulated_batches_ > 0) {
    double cell_volume = grid_->pitch() * grid_->pitch() * grid_->pitch();
    return it->second / n_accumulated_batches_ / cell_volume;
  }
  return 0.0;
}

//------------------------------------------------------------------------------

int FluxMesh::position_to_index(const std::array<double, 3>& position) const
{
  // 重要：必须与 FissionMatrix::position_to_index 使用完全相同的算法！
  // 否则同一个位置会被映射到不同的单元索引

  const auto& origin = grid_->origin();
  const auto& upper_bound = grid_->upper_bound();
  const auto& shape = grid_->shape();
  double inv_pitch = grid_->inv_pitch();

  int indices[3];
  for (int axis = 0; axis < 3; ++axis) {
    double coord = position[axis];
    // 边界检查 - 与 FissionMatrix 完全一致
    if (coord < origin[axis] || coord >= upper_bound[axis]) {
      return -1;
    }
    int idx = static_cast<int>((coord - origin[axis]) * inv_pitch);
    // 边界clamp - 与 FissionMatrix 完全一致
    if (idx < 0) {
      idx = 0;
    } else if (idx >= shape[axis]) {
      idx = shape[axis] - 1;
    }
    indices[axis] = idx;
  }

  // 线性索引计算 - 与 FissionMatrix 完全一致 (XYZ顺序)
  return (indices[0] * shape[1] + indices[1]) * shape[2] + indices[2];
}

//------------------------------------------------------------------------------

std::array<int, 3> FluxMesh::index_to_grid(int index) const
{
  // 逆向计算：index = (ix * shape[1] + iy) * shape[2] + iz
  const auto& shape = grid_->shape();
  std::array<int, 3> grid_idx;
  grid_idx[2] = index % shape[2]; // iz
  int remainder = index / shape[2];
  grid_idx[1] = remainder % shape[1]; // iy
  grid_idx[0] = remainder / shape[1]; // ix
  return grid_idx;
}

//------------------------------------------------------------------------------

int FluxMesh::grid_to_index(int ix, int iy, int iz) const
{
  // 线性索引计算 - 与 FissionMatrix 完全一致 (XYZ顺序)
  // index = (ix * shape[1] + iy) * shape[2] + iz
  const auto& shape = grid_->shape();
  return (ix * shape[1] + iy) * shape[2] + iz;
}

//------------------------------------------------------------------------------

void FluxMesh::write_hdf5(const std::string& filename, int n_batches) const
{
  // 只有主进程写入文件
#ifdef OPENMC_MPI
  if (mpi::rank != 0)
    return;
#endif

  if (n_accumulated_batches_ == 0) {
    std::cout << "FluxMesh: 无批次可写入,跳过文件输出." << std::endl;
    return;
  }

  // 创建 HDF5 文件
  hid_t file_id = file_open(filename, 'w');

  // 获取网格参数
  const auto& shape = grid_->shape();
  const auto& origin = grid_->origin();
  const auto& upper_bound = grid_->upper_bound();
  double pitch = grid_->pitch();
  double cell_volume = pitch * pitch * pitch;
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

  int n_groups = n_groups_;
  write_dataset(file_id, "n_groups", n_groups);
  if (!energy_edges_.empty()) {
    write_dataset(file_id, "energy_edges", energy_edges_);
  }

  // 准备稀疏数据数组
  std::vector<int> cell_indices;
  std::vector<double> flux_mean;
  std::vector<double> flux_std;
  std::vector<double> flux_group_mean;
  std::vector<double> flux_group_std;

  cell_indices.reserve(flux_sum_.size());
  flux_mean.reserve(flux_sum_.size());
  flux_std.reserve(flux_sum_.size());
  flux_group_mean.reserve(static_cast<size_t>(n_groups_) * flux_sum_.size());
  flux_group_std.reserve(static_cast<size_t>(n_groups_) * flux_sum_.size());

  // 计算每个非零单元的均值和标准差
  for (const auto& [cell_index, sum] : flux_sum_) {
    double mean = sum / n_accumulated_batches_ / cell_volume;
    double sum_sq = flux_sum_sq_.at(cell_index);
    double variance = 0.0;

    // 计算标准差: σ = sqrt(E[X²] - E[X]²)
    if (n_accumulated_batches_ > 1) {
      double mean_sq =
        sum_sq / n_accumulated_batches_ / (cell_volume * cell_volume);
      variance = mean_sq - mean * mean;
      // 防止数值误差导致负方差
      variance = std::max(variance, 0.0);
    }
    double std_dev = std::sqrt(variance);

    cell_indices.push_back(cell_index);
    flux_mean.push_back(mean);
    flux_std.push_back(std_dev);

    if (n_groups_ > 0) {
      std::vector<double> zero_vector = make_zero_group_vector();
      std::vector<double> zero_vector_sq = make_zero_group_vector();
      const auto& group_sum = flux_group_sum_.count(cell_index)
                                ? flux_group_sum_.at(cell_index)
                                : zero_vector;
      const auto& group_sum_sq = flux_group_sum_sq_.count(cell_index)
                                   ? flux_group_sum_sq_.at(cell_index)
                                   : zero_vector_sq;

      for (int g = 0; g < n_groups_; ++g) {
        double group_mean =
          group_sum[g] / n_accumulated_batches_ / cell_volume;
        double group_std = 0.0;
        if (n_accumulated_batches_ > 1) {
          double mean_sq = group_sum_sq[g] / n_accumulated_batches_ /
                           (cell_volume * cell_volume);
          double variance = mean_sq - group_mean * group_mean;
          group_std = std::sqrt(std::max(variance, 0.0));
        }
        flux_group_mean.push_back(group_mean);
        flux_group_std.push_back(group_std);
      }
    }
  }

  // 写入稀疏数据
  write_dataset(file_id, "cell_indices", cell_indices);
  write_dataset(file_id, "flux_mean", flux_mean);
  write_dataset(file_id, "flux_std", flux_std);
  if (n_groups_ > 0) {
    write_dataset(file_id, "flux_group_mean", flux_group_mean);
    write_dataset(file_id, "flux_group_std", flux_group_std);
  }

  // 创建密集格式的通量数组(用于可视化)
  std::vector<double> flux_dense(n_cells, 0.0);
  std::vector<double> flux_error_dense(n_cells, 0.0);
  std::vector<double> flux_group_mean_dense(
    static_cast<size_t>(n_groups_) * n_cells, 0.0);
  std::vector<double> flux_group_std_dense(
    static_cast<size_t>(n_groups_) * n_cells, 0.0);

  for (size_t i = 0; i < cell_indices.size(); ++i) {
    int idx = cell_indices[i];
    if (idx >= 0 && idx < static_cast<int>(n_cells)) {
      flux_dense[idx] = flux_mean[i];
      flux_error_dense[idx] = flux_std[i];

      if (n_groups_ > 0) {
        for (int g = 0; g < n_groups_; ++g) {
          size_t sparse_index = static_cast<size_t>(i) * n_groups_ + g;
          size_t dense_index = static_cast<size_t>(g) * n_cells + idx;
          if (sparse_index < flux_group_mean.size() &&
              dense_index < flux_group_mean_dense.size()) {
            flux_group_mean_dense[dense_index] = flux_group_mean[sparse_index];
            flux_group_std_dense[dense_index] = flux_group_std[sparse_index];
          }
        }
      }
    }
  }

  // 写入密集格式数据
  write_dataset(file_id, "flux_mean_dense", flux_dense);
  write_dataset(file_id, "flux_std_dense", flux_error_dense);
  if (n_groups_ > 0) {
    write_dataset(file_id, "flux_group_mean_dense", flux_group_mean_dense);
    write_dataset(file_id, "flux_group_std_dense", flux_group_std_dense);
  }

  // 关闭文件
  file_close(file_id);
}

std::vector<double> FluxMesh::make_zero_group_vector() const
{
  return std::vector<double>(static_cast<size_t>(n_groups_), 0.0);
}

int FluxMesh::determine_group(double energy_eV, int mg_group) const
{
  if (n_groups_ <= 1) {
    return 0;
  }

  if (mg_group >= 0) {
    if (mg_group >= n_groups_) {
      return n_groups_ - 1;
    }
    return mg_group;
  }

  if (!energy_edges_.empty() && energy_eV >= 0.0) {
    if (energy_eV < energy_edges_.front()) {
      return 0;
    }
    if (energy_eV >= energy_edges_.back()) {
      return n_groups_ - 1;
    }

    auto it =
      std::upper_bound(energy_edges_.begin(), energy_edges_.end(), energy_eV);
    int idx = static_cast<int>(std::distance(energy_edges_.begin(), it)) - 1;
    if (idx < 0) {
      idx = 0;
    }
    if (idx >= n_groups_) {
      idx = n_groups_ - 1;
    }
    return idx;
  }

  return -1;
}

} // namespace openmc
