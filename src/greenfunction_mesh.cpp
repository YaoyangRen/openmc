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
  upper_bound_ = grid_->upper_bound();

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
  // 重要：必须与 FissionMatrix::position_to_index 使用完全相同的算法！
  // 否则同一个位置会被映射到不同的单元索引，导致伴随源覆盖不足

  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();
  const auto& upper_bound = upper_bound_;

  int indices[3];
  double coords[3] = {r.x, r.y, r.z};

  for (int axis = 0; axis < 3; ++axis) {
    double coord = coords[axis];
    // 边界检查 - 与 FissionMatrix 完全一致
    if (coord < origin[axis] || coord >= upper_bound[axis]) {
      return -1;
    }
    int idx = static_cast<int>((coord - origin[axis]) * inv_pitch_);
    // 边界clamp - 与 FissionMatrix 完全一致
    if (idx < 0) {
      idx = 0;
    } else if (idx >= shape[axis]) {
      idx = shape[axis] - 1;
    }
    indices[axis] = idx;
  }

  // 线性索引计算 - 与 FissionMatrix 完全一致
  return (indices[0] * shape[1] + indices[1]) * shape[2] + indices[2];
}

void GreenFunctionMesh::record_source_birth(
  const Position& r, int64_t source_particle_id, double energy, int mg_group)
{
  // 验证源粒子ID
  if (source_particle_id < 0) {
    return;
  }

  // 计算源单元索引
  int i_source_cell = position_to_cell_index(r);

  if (i_source_cell < 0) {
    // 源粒子在网格边界外，静默忽略
    return;
  }

  // 确定源能群
  int g_source = determine_group(energy, mg_group);
  // source_state = source_cell * n_source_groups + g_source
  int source_state = i_source_cell * n_groups_ + g_source;

  // 记录源粒子到源状态的映射
  particle_to_source_state_[source_particle_id] = source_state;

  // 增加该源单元的计数（仍然按 cell 计数，用于 source_counts_ 输出）
  source_counts_[i_source_cell]++;
}

void GreenFunctionMesh::accumulate(const Position& r, double contribution,
  int64_t source_particle_id, int family, double energy_eV, int mg_group)
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

  // 查找源状态 (source_state = source_cell * n_source_groups + g_source)
  auto it = particle_to_source_state_.find(source_particle_id);
  if (it == particle_to_source_state_.end()) {
    // 这个粒子没有记录源位置，忽略
    return;
  }

  int i_source = it->second; // source_state

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

  // Clamp family index: 0=prompt, 1..8=delayed
  if (family < 0 || family >= N_FAMILIES)
    family = 0;
  size_t flat_idx = static_cast<size_t>(family) * n_groups_ + group;

  // 稀疏存储：只在有贡献时才创建条目
  {
    std::lock_guard<std::mutex> lock(data_mutex_);

    auto& response_map = current_batch_transfer_data_sparse_[i_source];
    auto& group_vector = response_map[j_response];
    if (group_vector.empty()) {
      group_vector = make_zero_group_vector();
    }
    group_vector[flat_idx] += contribution;

    auto& cumulative_vector = cumulative_data_sparse_[j_response];
    if (cumulative_vector.empty()) {
      cumulative_vector = make_zero_group_vector();
    }
    cumulative_vector[flat_idx] += contribution;
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
        for (size_t g = 0; g < values.size(); ++g) {
          target_vector[g] += values[g];
        }
      }
    }
  }

  // 开始新batch
  current_batch_id_ = batch_id;
  current_batch_transfer_data_sparse_.clear();
  // 清空源粒子→源状态映射，避免内存持续增长
  particle_to_source_state_.clear();
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

  // 计算有源的状态数和总源粒子数
  // transfer_functions_sparse_ 外层 key 现在是 source_state
  int n_source_states = transfer_functions_sparse_.size();
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
                                     (n_source_states * spatial_size_));
  size_t values_per_entry = static_cast<size_t>(N_FAMILIES) * n_groups_;
  size_t sparse_memory_mb =
    (total_nonzero_entries * (sizeof(int) + values_per_entry * sizeof(double)) +
      n_source_states * 64) /
    (1024 * 1024);
  size_t dense_memory_mb =
    (n_source_states * spatial_size_ * sizeof(double)) / (1024 * 1024);

  // 简洁输出传递函数信息
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "传递函数 T(source_state->r): " << n_source_states
            << " 个源状态 (" << n_groups_ << " 群×cells), "
            << total_source_particles << " 个源粒子, " << total_contributions_
            << " 次贡献" << std::endl;
  std::cout << "  稀疏存储: " << total_nonzero_entries << " 个非零条目 ("
            << std::fixed << std::setprecision(1) << (100.0 - sparsity)
            << "% 密度)" << std::endl;
  if (dense_memory_mb > 0) {
    std::cout << "  内存节省: " << dense_memory_mb << " MB -> "
              << sparse_memory_mb << " MB (节省 " << std::setprecision(0)
              << (100.0 * (1.0 - static_cast<double>(sparse_memory_mb) /
                                   dense_memory_mb))
              << "%)" << std::endl;
  }
  std::cout << "  输出: " << filename << std::endl;

  // 创建HDF5文件
  hid_t file_id = file_open(filename, 'w');

  // 写入文件头部信息
  write_attribute(file_id, "filetype", "transfer_function_mesh_sparse");
  write_attribute(file_id, "version", "5.0"); // 新版本：source_state 扩展
  write_attribute(file_id, "description",
    "Transfer function T(source_state->r_response): "
    "source_state = source_cell * n_source_groups + g_source");
  write_attribute(file_id, "pitch", grid_->pitch());
  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());
  std::array<double, 3> grid_pitch {
    grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(file_id, "grid_shape", grid_->shape());
  write_dataset(file_id, "grid_lower_left", grid_->origin());
  write_dataset(file_id, "grid_upper_right", grid_->upper_bound());
  write_dataset(file_id, "grid_pitch", grid_pitch);
  write_attribute(file_id, "n_source_states", n_source_states);
  write_attribute(file_id, "n_source_groups", n_groups_);
  write_attribute(file_id, "total_source_particles", total_source_particles);
  write_attribute(file_id, "total_nonzero_entries",
    static_cast<int64_t>(total_nonzero_entries));
  write_attribute(file_id, "sparsity_percent", sparsity);
  write_dataset(file_id, "n_groups", n_groups_);
  write_dataset(file_id, "n_families", N_FAMILIES);
  write_attribute(file_id, "family_order", "prompt,delayed_1,...,delayed_8");
  write_attribute(file_id, "storage_order", "family_major: [family][group]");
  if (!energy_edges_.empty()) {
    write_dataset(file_id, "energy_edges", energy_edges_);
  }

  // 写入累积的传递函数（稀疏格式）
  vector<int> cumulative_indices;
  vector<double> cumulative_values;
  cumulative_indices.reserve(cumulative_data_sparse_.size());
  cumulative_values.reserve(cumulative_data_sparse_.size() *
                            static_cast<size_t>(N_FAMILIES) * n_groups_);

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

  // 收集有源的状态索引列表（排序）
  // source_state = source_cell * n_source_groups + g_source
  vector<int> source_state_indices;
  source_state_indices.reserve(n_source_states);
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    source_state_indices.push_back(i_source);
  }
  std::sort(source_state_indices.begin(), source_state_indices.end());
  write_dataset(file_id, "source_state_indices", source_state_indices);

  // 同时写出每个源状态对应的 cell 和 group（方便后处理解码）
  {
    vector<int> ss_cells, ss_groups;
    ss_cells.reserve(n_source_states);
    ss_groups.reserve(n_source_states);
    for (int ss : source_state_indices) {
      ss_cells.push_back(ss / n_groups_);
      ss_groups.push_back(ss % n_groups_);
    }
    write_dataset(file_id, "source_state_cell", ss_cells);
    write_dataset(file_id, "source_state_group", ss_groups);
  }

  // 仍然输出 source_cell_indices（向后兼容），通过唯一化 source_state_cell 得到
  {
    vector<int> cell_ids;
    for (const auto& count : source_counts_) {
      if (count > 0) {
        cell_ids.push_back(static_cast<int>(&count - source_counts_.data()));
      }
    }
    if (!cell_ids.empty())
      write_dataset(file_id, "source_cell_indices", cell_ids);
  }

  // 写入每个源单元的粒子计数
  write_dataset(file_id, "source_counts_per_cell", source_counts_);

  // 为每个源状态创建一个组（稀疏格式）
  hid_t transfer_group = H5Gcreate(
    file_id, "transfer_functions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // 写入每个源状态的传递函数 T(source_state -> r)
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    std::string state_name = "source_state_" + std::to_string(i_source);
    hid_t cell_group = H5Gcreate(transfer_group, state_name.c_str(),
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 收集该源状态的非零响应
    vector<int> response_indices;
    vector<double> response_values;
    response_indices.reserve(response_map.size());
    response_values.reserve(
      response_map.size() * static_cast<size_t>(N_FAMILIES) * n_groups_);

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
    write_attribute(cell_group, "source_cell", i_source / n_groups_);
    write_attribute(cell_group, "source_group", i_source % n_groups_);

    H5Gclose(cell_group);
  }

  H5Gclose(transfer_group);
  file_close(file_id);

  // 清理数据
  transfer_functions_sparse_.clear();
  current_batch_transfer_data_sparse_.clear();
  cumulative_data_sparse_.clear();
  source_counts_.clear();
  particle_to_source_state_.clear();
}

std::vector<double> GreenFunctionMesh::make_zero_group_vector() const
{
  return std::vector<double>(static_cast<size_t>(N_FAMILIES) * n_groups_, 0.0);
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
