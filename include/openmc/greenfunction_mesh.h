#ifndef OPENMC_GREENFUNCTION_MESH_H
#define OPENMC_GREENFUNCTION_MESH_H

#include "hdf5.h"
#include "openmc/array.h"
#include "openmc/mesh_init.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

//! \class GreenFunctionMesh
//! \brief 传递函数（Transfer Function）网格：T(r_source -> r_response)
//!
//! 根据传递函数理论：
//!   T(r_s -> r) = ∫∫ ν̄Σf(r,E')Φ(r,Ω',E'|S(r_s)) dΩ'dE'
//!
//! 物理意义：
//!   从源位置 r_s 出发的中子在响应位置 r 处产生的
//!   平均裂变中子数的期望值
//!
//! 数据结构：
//!   - 按源单元索引存储：map<i_source, T(i -> r)>
//!   - 便于与伴随源 S†(r_s) 进行空间卷积
//!   - 共轭通量：Φ†(r) = Σ_i T(i -> r) × S†(i)
//!
//! 应用：
//!   - 计算共轭通量（伴随通量）
//!   - 重要性函数求解
//!   - 扰动分析
//!
class GreenFunctionMesh {
public:
  //! prompt(0) + delayed_1..8 共 9 个 family
  static constexpr int N_FAMILIES = 9;

  // 构造函数，初始化传递函数网格
  explicit GreenFunctionMesh(std::shared_ptr<SharedMeshGrid> grid,
    int max_batches, std::vector<double> energy_edges = {});

  // 记录源粒子的出生位置（建立源粒子ID到源单元的映射）
  void record_source_birth(const Position& r, int64_t source_particle_id);

  // 为特定源粒子累积传递函数贡献
  // contribution: 期望裂变中子数 nu_t = (w/k_eff) × w_ufs × (ν̄Σf/Σt)
  // family: 0=prompt, 1..8=delayed_group (来自 p.delayed_group())
  void accumulate(const Position& r, double contribution,
    int64_t source_particle_id, int family = 0, double energy_eV = -1.0,
    int mg_group = -1);

  // 开始新batch
  void start_new_batch(int batch_id);

  // 写入文件（可指定文件名）
  void finalize_greenfunction_mesh(const int batch_id,
    const std::string& filename = "transfer_function_data.h5");

  // 获取特定源单元的传递函数数据 T(i_source -> r)
  const vector<double>& get_source_cell_data(int source_cell_index) const;

  // 获取源单元统计信息
  const vector<int>& get_source_counts() const { return source_counts_; }
  const vector<int>& get_source_cell_indices() const;

  // 网格尺寸信息
  const std::array<int, 3>& shape() const { return grid_->shape(); }
  const std::array<double, 3>& origin() const { return grid_->origin(); }
  double pitch() const { return grid_->pitch(); }

  // 计算空间位置对应的单元索引
  int position_to_cell_index(const Position& r) const;

private:
  // 核心数据：按源单元索引存储的传递函数（稀疏存储）
  // 外层 Key: 源单元索引 i_source (0 到 nx*ny*nz-1)
  // 内层 Key: 响应单元索引 j_response (0 到 nx*ny*nz-1)
  // Value: T(i_source -> j_response) 传递函数值
  std::unordered_map<int, std::unordered_map<int, std::vector<double>>>
    transfer_functions_sparse_;

  // 当前batch中每个源单元的传递函数数据（稀疏）
  std::unordered_map<int, std::unordered_map<int, std::vector<double>>>
    current_batch_transfer_data_sparse_;

  // 每个源单元产生的源粒子计数 [nx*ny*nz]
  vector<int> source_counts_;

  // 源粒子ID到源单元索引的映射
  std::unordered_map<int64_t, int> particle_to_source_cell_;

  // 累积所有源单元的传递函数（总的传递函数，稀疏存储）
  std::unordered_map<int, std::vector<double>> cumulative_data_sparse_;

  // 统一网格配置
  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_; // 上边界（用于输出）
  double inv_pitch_;                  // 网格单元边长的倒数
  int current_batch_id_;              // 当前处理的batch ID
  int max_batches_;                   // 最大batch数量
  size_t spatial_size_;               // 统计信息
  std::vector<double> energy_edges_;
  int n_groups_ {1};
  std::atomic<uint64_t> dropped_contributions_ {0};
  std::atomic<uint64_t> total_contributions_ {0};

  // 线程同步
  mutable std::mutex data_mutex_; // 保护current_batch_particle_data_的访问

  // 分群辅助函数
  std::vector<double> make_zero_group_vector() const;
  int determine_group(double energy_eV, int mg_group) const;
};

} // namespace openmc

#endif // OPENMC_GREENFUNCTION_MESH_H