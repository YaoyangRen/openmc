#include "openmc/adjoint_flux.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"

namespace openmc {

void AdjointFlux::compute_from_files(const std::string& transfer_function_file,
  const std::string& fission_matrix_file, const std::string& output_file)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "ADJOINT FLUX COMPUTATION" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  // 1. 读取传递函数数据
  std::cout << "\nReading transfer function data from: "
            << transfer_function_file << std::endl;

  // 检查传递函数文件是否存在
  if (!file_exists(transfer_function_file)) {
    fatal_error(
      "Transfer function file not found: " + transfer_function_file +
      "\nPlease ensure the transfer function has been computed and saved.");
  }

  hid_t tf_file = file_open(transfer_function_file, 'r');

  // 读取网格参数
  read_dataset(tf_file, "origin", origin_);
  read_dataset(tf_file, "shape", shape_);
  read_attribute(tf_file, "pitch", pitch_);

  n_cells_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  std::cout << "  Grid: " << shape_[0] << " x " << shape_[1] << " x "
            << shape_[2] << " = " << n_cells_ << " cells" << std::endl;
  std::cout << "  Pitch: " << pitch_ << " cm" << std::endl;

  // 读取源单元索引列表
  vector<int> source_cell_indices;
  read_dataset(tf_file, "source_cell_indices", source_cell_indices);

  std::cout << "  Source cells: " << source_cell_indices.size() << std::endl;

  // 读取传递函数（稀疏格式）
  std::unordered_map<int, std::unordered_map<int, double>> transfer_functions;

  hid_t tf_group = H5Gopen(tf_file, "transfer_functions", H5P_DEFAULT);

  size_t total_tf_entries = 0;
  for (int i_source : source_cell_indices) {
    std::string cell_name = "source_cell_" + std::to_string(i_source);

    if (object_exists(tf_group, cell_name.c_str())) {
      hid_t cell_group = H5Gopen(tf_group, cell_name.c_str(), H5P_DEFAULT);

      vector<int> indices;
      vector<double> values;
      read_dataset(cell_group, "indices", indices);
      read_dataset(cell_group, "values", values);

      // 构建稀疏 map
      for (size_t k = 0; k < indices.size(); ++k) {
        transfer_functions[i_source][indices[k]] = values[k];
      }

      total_tf_entries += indices.size();
      H5Gclose(cell_group);
    }
  }

  H5Gclose(tf_group);
  file_close(tf_file);

  std::cout << "  Transfer function entries: " << total_tf_entries << std::endl;

  // 2. 读取伴随源数据
  std::cout << "\nReading adjoint source from: " << fission_matrix_file
            << std::endl;

  // 检查裂变矩阵文件是否存在
  if (!file_exists(fission_matrix_file)) {
    fatal_error(
      "Fission matrix file not found: " + fission_matrix_file +
      "\nPlease ensure the fission matrix has been computed and saved.");
  }

  hid_t fm_file = file_open(fission_matrix_file, 'r');

  // 检查是否存在伴随源
  if (!object_exists(fm_file, "adjoint_source")) {
    file_close(fm_file);
    fatal_error("Adjoint source not found in fission matrix file. Run "
                "compute_adjoint_source first.");
  }

  vector<double> adjoint_source;
  read_dataset(fm_file, "adjoint_source", adjoint_source);

  double k_adjoint = 0.0;
  if (attribute_exists(fm_file, "k_adjoint")) {
    read_attribute(fm_file, "k_adjoint", k_adjoint);
  }

  file_close(fm_file);

  std::cout << "  Adjoint source cells: " << adjoint_source.size() << std::endl;
  std::cout << "  k_adjoint: " << std::fixed << std::setprecision(6)
            << k_adjoint << std::endl;

  // 验证网格一致性
  if (adjoint_source.size() != n_cells_) {
    fatal_error("Grid mismatch: transfer function has " +
                std::to_string(n_cells_) + " cells, but adjoint source has " +
                std::to_string(adjoint_source.size()) + " cells");
  }

  // 3. 计算共轭通量
  std::cout << "\nComputing adjoint flux: Φ†(r) = Σ_i T(i→r) × S†(i)"
            << std::endl;

  compute_from_memory(
    transfer_functions, adjoint_source, shape_, origin_, pitch_);

  // 4. 写入输出文件
  write_to_file(output_file);

  std::cout << std::string(70, '=') << std::endl;
}

void AdjointFlux::compute_from_memory(
  const std::unordered_map<int, std::unordered_map<int, double>>&
    transfer_functions,
  const vector<double>& adjoint_source, const std::array<int, 3>& shape,
  const std::array<double, 3>& origin, double pitch)
{
  shape_ = shape;
  origin_ = origin;
  pitch_ = pitch;
  n_cells_ = static_cast<size_t>(shape[0]) * shape[1] * shape[2];

  // 清空之前的结果
  adjoint_flux_sparse_.clear();

  // 计算 Φ†(r) = Σ_i T(i -> r) × S†(i)
  for (const auto& [i_source, response_map] : transfer_functions) {
    // 获取源单元的伴随源值（重要性）
    if (i_source < 0 || i_source >= static_cast<int>(adjoint_source.size())) {
      continue;
    }

    double importance = adjoint_source[i_source];

    if (importance <= 0.0) {
      continue; // 跳过零重要性的源
    }

    // 对该源单元的所有响应位置进行累加
    for (const auto& [j_response, T_value] : response_map) {
      // Φ†(j) += T(i -> j) × S†(i)
      adjoint_flux_sparse_[j_response] += T_value * importance;
    }
  }

  // 计算统计信息
  nonzero_cells_ = adjoint_flux_sparse_.size();
  max_flux_ = 0.0;
  total_flux_ = 0.0;

  for (const auto& [j, flux] : adjoint_flux_sparse_) {
    max_flux_ = std::max(max_flux_, flux);
    total_flux_ += flux;
  }

  // 输出统计信息
  double density = 100.0 * nonzero_cells_ / n_cells_;
  std::cout << "\nAdjoint Flux Statistics:" << std::endl;
  std::cout << "  Nonzero cells: " << nonzero_cells_ << " / " << n_cells_
            << " (" << std::fixed << std::setprecision(2) << density << "%)"
            << std::endl;
  std::cout << "  Max value: " << std::scientific << std::setprecision(6)
            << max_flux_ << std::endl;
  std::cout << "  Total flux: " << total_flux_ << std::endl;
  std::cout << "  Mean (nonzero): "
            << (nonzero_cells_ > 0 ? total_flux_ / nonzero_cells_ : 0.0)
            << std::endl;
}

vector<double> AdjointFlux::get_adjoint_flux_dense() const
{
  vector<double> dense(n_cells_, 0.0);

  for (const auto& [j, flux] : adjoint_flux_sparse_) {
    if (j >= 0 && j < static_cast<int>(n_cells_)) {
      dense[j] = flux;
    }
  }

  return dense;
}

double AdjointFlux::get_max_value() const
{
  return max_flux_;
}

double AdjointFlux::get_total_flux() const
{
  return total_flux_;
}

void AdjointFlux::write_to_file(const std::string& filename)
{
  std::cout << "\nWriting adjoint flux to: " << filename << std::endl;

  hid_t file_id = file_open(filename, 'w');

  // 写入文件属性
  write_attribute(file_id, "filetype", "adjoint_flux");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "description",
    "Adjoint flux computed from transfer function and adjoint source");
  write_attribute(file_id, "storage_format", "sparse");
  write_attribute(file_id, "pitch", pitch_);

  // 写入网格信息
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells_));

  // 写入统计信息
  write_attribute(file_id, "nonzero_cells", static_cast<int>(nonzero_cells_));
  write_attribute(file_id, "max_flux", max_flux_);
  write_attribute(file_id, "total_flux", total_flux_);
  write_attribute(
    file_id, "density_percent", 100.0 * nonzero_cells_ / n_cells_);

  // 写入稀疏格式的共轭通量数据
  vector<int> indices;
  vector<double> values;

  indices.reserve(adjoint_flux_sparse_.size());
  values.reserve(adjoint_flux_sparse_.size());

  for (const auto& [j, flux] : adjoint_flux_sparse_) {
    indices.push_back(j);
    values.push_back(flux);
  }

  // 按索引排序
  vector<size_t> sort_indices(indices.size());
  std::iota(sort_indices.begin(), sort_indices.end(), 0);
  std::sort(sort_indices.begin(), sort_indices.end(),
    [&indices](size_t a, size_t b) { return indices[a] < indices[b]; });

  vector<int> sorted_indices(indices.size());
  vector<double> sorted_values(values.size());
  for (size_t i = 0; i < sort_indices.size(); ++i) {
    sorted_indices[i] = indices[sort_indices[i]];
    sorted_values[i] = values[sort_indices[i]];
  }

  write_dataset(file_id, "cell_indices", sorted_indices);
  write_dataset(file_id, "adjoint_flux_values", sorted_values);

  // 也写入稠密格式（可选，用于可视化）
  auto dense_flux = get_adjoint_flux_dense();
  write_dataset(file_id, "adjoint_flux_dense", dense_flux);

  file_close(file_id);

  std::cout << "  Sparse entries written: " << adjoint_flux_sparse_.size()
            << std::endl;
  std::cout << "  File size estimate: "
            << (adjoint_flux_sparse_.size() * 12 + n_cells_ * 8) / (1024 * 1024)
            << " MB" << std::endl;
}

} // namespace openmc
