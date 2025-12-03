#include "openmc/greenfunction_mesh.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/universe.h"
#include "openmc/vector.h"

namespace openmc {

GreenFunctionMesh::GreenFunctionMesh(std::shared_ptr<SharedMeshGrid> grid,
  int max_batches, std::vector<double> energy_edges)
  : grid_(grid), inv_pitch_(1.0 / grid->pitch()), current_batch_id_(-1),
    max_batches_(max_batches), energy_edges_(std::move(energy_edges))
{
  // 计算上边界（用于输出）
  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();
  double pitch = grid_->pitch();
  upper_bound_ = {origin[0] + shape[0] * pitch, origin[1] + shape[1] * pitch,
    origin[2] + shape[2] * pitch};

  spatial_size_ = grid_->n_cells();

  if (energy_edges_.size() >= 2) {
    if (!std::is_sorted(energy_edges_.begin(), energy_edges_.end())) {
      fatal_error(
        "GreenFunctionMesh: energy group edges must be sorted ascending.");
    }
    n_groups_ = static_cast<int>(energy_edges_.size()) - 1;
  } else {
    energy_edges_.clear();
    n_groups_ = 1;
  }

  // 初始化信息将在finalize时输出

  // 初始化源计数
  source_counts_.resize(spatial_size_, 0);

  // 稀疏存储：transfer_functions_sparse_ 和 cumulative_data_sparse_ 按需增长
}

int GreenFunctionMesh::position_to_cell_index(const Position& r) const
{
  constexpr double eps = 1.0e-10;

  int ix =
    static_cast<int>(std::floor((r.x - grid_->origin()[0] + eps) * inv_pitch_));
  int iy =
    static_cast<int>(std::floor((r.y - grid_->origin()[1] + eps) * inv_pitch_));
  int iz =
    static_cast<int>(std::floor((r.z - grid_->origin()[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix < 0 || ix >= grid_->shape()[0] || iy < 0 || iy >= grid_->shape()[1] ||
      iz < 0 || iz >= grid_->shape()[2]) {
    return -1; // 超出边界
  }

  // 计算线性索引
  return ix * grid_->shape()[1] * grid_->shape()[2] + iy * grid_->shape()[2] +
         iz;
}

void GreenFunctionMesh::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  // 验证源粒子ID
  if (source_particle_id < 0) {
    return;
  }

  // 计算源单元索引
  int i_source = position_to_cell_index(r);

  if (i_source < 0) {
    // 源粒子在网格边界外，静默忽略
    return;
  }

  // 记录源粒子到源单元的映射
  particle_to_source_cell_[source_particle_id] = i_source;

  // 增加该源单元的计数
  source_counts_[i_source]++;

  // 稀疏存储：不需要预分配，map 会按需增长
  // transfer_functions_sparse_[i_source] 和
  // current_batch_transfer_data_sparse_[i_source] 会在 accumulate()
  // 中首次使用时自动创建
}

void GreenFunctionMesh::accumulate(const Position& r, double contribution,
  int64_t source_particle_id, double energy_eV, int mg_group)
{
  // 累积传递函数 T(r_source -> r_response):
  // contribution = nu_t = (w/k_eff) × w_ufs × (ν̄Σf/Σt)
  // 表示源位置在响应位置 r 处产生的期望裂变中子数

  // 验证源粒子ID
  if (source_particle_id < 0) {
    return;
  }

  // 验证贡献值（稀疏存储：只存储非零值）
  if (!std::isfinite(contribution) || contribution <= 0.0) {
    return;
  }

  // 查找源单元
  auto it = particle_to_source_cell_.find(source_particle_id);
  if (it == particle_to_source_cell_.end()) {
    // 这个粒子没有记录源位置，忽略
    return;
  }

  int i_source = it->second;

  // 计算响应位置的单元索引
  int j_response = position_to_cell_index(r);

  if (j_response < 0) {
    dropped_contributions_++;
    return;
  }

  int group = determine_group(energy_eV, mg_group);
  if (group < 0) {
    fatal_error(
      "GreenFunctionMesh: unable to determine energy group for contribution.");
  }

  // 稀疏存储：只在有贡献时才创建条目
  {
    std::lock_guard<std::mutex> lock(data_mutex_);

    auto& response_map = current_batch_transfer_data_sparse_[i_source];
    auto& group_vector = response_map[j_response];
    if (group_vector.empty()) {
      group_vector = make_zero_group_vector();
    }
    group_vector[group] += contribution;

    auto& cumulative_vector = cumulative_data_sparse_[j_response];
    if (cumulative_vector.empty()) {
      cumulative_vector = make_zero_group_vector();
    }
    cumulative_vector[group] += contribution;
  }

  total_contributions_++;
}
void GreenFunctionMesh::start_new_batch(int batch_id)
{
  // 保存上一个batch的数据（稀疏）
  if (current_batch_id_ >= 0) {
    for (const auto& [i_source, response_map] :
      current_batch_transfer_data_sparse_) {
      auto& target_map = transfer_functions_sparse_[i_source];
      for (const auto& [j_response, values] : response_map) {
        auto& target_vector = target_map[j_response];
        if (target_vector.empty()) {
          target_vector = make_zero_group_vector();
        }
        for (int g = 0; g < n_groups_; ++g) {
          target_vector[g] += values[g];
        }
      }
    }
  }

  // 开始新batch
  current_batch_id_ = batch_id;
  current_batch_transfer_data_sparse_.clear();
}

const vector<double>& GreenFunctionMesh::get_source_cell_data(
  int source_cell_index) const
{
  // 从稀疏存储构建稠密向量
  static vector<double> dense_data;

  auto it = transfer_functions_sparse_.find(source_cell_index);
  if (it != transfer_functions_sparse_.end()) {
    dense_data.clear();
    dense_data.resize(spatial_size_, 0.0);

    // 填充非零值（按群求和得到总贡献）
    for (const auto& [j_response, values] : it->second) {
      if (j_response >= 0 && j_response < static_cast<int>(spatial_size_)) {
        double total = std::accumulate(values.begin(), values.end(), 0.0);
        dense_data[j_response] = total;
      }
    }
    return dense_data;
  }

  // 返回空向量
  static vector<double> empty_data;
  return empty_data;
}

const vector<int>& GreenFunctionMesh::get_source_cell_indices() const
{
  static vector<int> indices;
  indices.clear();
  indices.reserve(transfer_functions_sparse_.size());

  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    indices.push_back(i_source);
  }

  std::sort(indices.begin(), indices.end());
  return indices;
}
void GreenFunctionMesh::finalize_greenfunction_mesh(
  const int batch_id, const std::string& filename)
{
  // 保存最后一个batch的数据
  start_new_batch(-1);

  if (transfer_functions_sparse_.empty()) {
    return;
  }

  // 计算有源的单元数和总源粒子数
  int n_source_cells = transfer_functions_sparse_.size();
  int total_source_particles = 0;
  for (const auto& count : source_counts_) {
    total_source_particles += count;
  }

  // 计算稀疏性统计
  size_t total_nonzero_entries = 0;
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    total_nonzero_entries += response_map.size();
  }

  double sparsity = 100.0 * (1.0 - static_cast<double>(total_nonzero_entries) /
                                     (n_source_cells * spatial_size_));
  size_t sparse_memory_mb =
    (total_nonzero_entries * (sizeof(int) + sizeof(double)) +
      n_source_cells * 64) /
    (1024 * 1024); // 估算
  size_t dense_memory_mb =
    (n_source_cells * spatial_size_ * sizeof(double)) / (1024 * 1024);

  // 简洁输出传递函数信信息
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "\nTransfer Function T(r_s->r): " << n_source_cells
            << " source cells, " << total_source_particles
            << " source particles, " << total_contributions_ << " contributions"
            << std::endl;
  std::cout << "  Sparse storage: " << total_nonzero_entries
            << " non-zero entries (" << std::fixed << std::setprecision(2)
            << (100.0 - sparsity) << "% density, " << sparsity << "% sparsity)"
            << std::endl;
  std::cout << "  Memory saved: " << dense_memory_mb << " MB (dense) -> "
            << sparse_memory_mb << " MB (sparse), " << std::fixed
            << std::setprecision(1)
            << (100.0 * (1.0 - static_cast<double>(sparse_memory_mb) /
                                 dense_memory_mb))
            << "% reduction" << std::endl;
  std::cout << "  Output: " << filename << std::endl;

  // 创建HDF5文件，使用指定的文件名
  hid_t file_id = file_open(filename, 'w');

  // 写入文件头部信息
  write_attribute(file_id, "filetype", "transfer_function_mesh_sparse");
  write_attribute(file_id, "version", "3.0");
  write_attribute(file_id, "description",
    "Transfer function T(r_source->r_response): Sparse storage format");
  write_attribute(file_id, "pitch", grid_->pitch());
  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());
  write_attribute(file_id, "n_source_cells", n_source_cells);
  write_attribute(file_id, "total_source_particles", total_source_particles);
  write_attribute(file_id, "total_nonzero_entries",
    static_cast<int64_t>(total_nonzero_entries));
  write_attribute(file_id, "sparsity_percent", sparsity);
  write_dataset(file_id, "n_groups", n_groups_);
  if (!energy_edges_.empty()) {
    write_dataset(file_id, "energy_edges", energy_edges_);
  }

  // 写入累积的传递函数（稀疏格式）
  vector<int> cumulative_indices;
  vector<double> cumulative_values;
  cumulative_indices.reserve(cumulative_data_sparse_.size());
  cumulative_values.reserve(
    cumulative_data_sparse_.size() * static_cast<size_t>(n_groups_));

  for (const auto& [j_response, values] : cumulative_data_sparse_) {
    cumulative_indices.push_back(j_response);
    cumulative_values.insert(
      cumulative_values.end(), values.begin(), values.end());
  }

  hid_t cumulative_group = H5Gcreate(file_id, "cumulative_transfer_function",
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  write_dataset(cumulative_group, "indices", cumulative_indices);
  write_dataset(cumulative_group, "values", cumulative_values);
  H5Gclose(cumulative_group);

  // 收集有源的单元索引列表（排序�?
  vector<int> source_cell_indices;
  source_cell_indices.reserve(n_source_cells);
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    source_cell_indices.push_back(i_source);
  }
  std::sort(source_cell_indices.begin(), source_cell_indices.end());
  write_dataset(file_id, "source_cell_indices", source_cell_indices);

  // 收集每个源单元的粒子计数（仅有源的单元）
  vector<int> source_counts_nonzero;
  source_counts_nonzero.reserve(n_source_cells);
  for (int i_source : source_cell_indices) {
    source_counts_nonzero.push_back(source_counts_[i_source]);
  }
  write_dataset(file_id, "source_counts", source_counts_nonzero);

  // 为每个源单元创建一个组（稀疏格式）
  hid_t transfer_group = H5Gcreate(
    file_id, "transfer_functions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // 写入每个源单元的传递函�?T(i_source -> r) - 稀疏格�?
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    std::string cell_name = "source_cell_" + std::to_string(i_source);
    hid_t cell_group = H5Gcreate(
      transfer_group, cell_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 收集该源单元的非零响�?
    vector<int> response_indices;
    vector<double> response_values;
    response_indices.reserve(response_map.size());
    response_values.reserve(response_map.size() * n_groups_);

    for (const auto& [j_response, values] : response_map) {
      response_indices.push_back(j_response);
      response_values.insert(
        response_values.end(), values.begin(), values.end());
    }

    write_dataset(cell_group, "indices", response_indices);
    write_dataset(cell_group, "values", response_values);
    write_attribute(
      cell_group, "n_nonzero", static_cast<int>(response_map.size()));
    write_attribute(cell_group, "n_groups", n_groups_);

    H5Gclose(cell_group);
  }

  H5Gclose(transfer_group);
  file_close(file_id);

  // 清理数据
  transfer_functions_sparse_.clear();
  current_batch_transfer_data_sparse_.clear();
  cumulative_data_sparse_.clear();
  source_counts_.clear();
  particle_to_source_cell_.clear();
}

std::vector<double> GreenFunctionMesh::make_zero_group_vector() const
{
  return std::vector<double>(static_cast<size_t>(n_groups_), 0.0);
}

int GreenFunctionMesh::determine_group(double energy_eV, int mg_group) const
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
