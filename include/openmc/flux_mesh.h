#ifndef OPENMC_FLUX_MESH_H
#define OPENMC_FLUX_MESH_H

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

// 前向声明
class SharedMeshGrid;

//==============================================================================
//! 通量网格类 - 用于统计空间通量分布
//! 使用与传递函数和裂变矩阵相同的空间网格
//==============================================================================
class FluxMesh {
public:
  //! 构造函数
  //! \param grid 共享的网格对象
  explicit FluxMesh(std::shared_ptr<SharedMeshGrid> grid,
    std::vector<double> energy_edges = {});

  //! 析构函数
  ~FluxMesh() = default;

  //! 累积单个粒子对通量的贡献
  //! \param position 粒子位置 [x, y, z]
  //! \param weight 粒子权重
  //! \param distance 径迹长度(用于径迹长度估计器)
  void accumulate(const std::array<double, 3>& position, double weight,
    double distance, double energy_eV = -1.0, int mg_group = -1);

  //! 批次结束处理 - 累积批次统计
  //! \param batch 当前批次号
  void end_batch(int batch);

  //! 模拟结束处理 - 输出结果到 HDF5 文件
  //! \param n_batches 总批次数
  void finalize(int n_batches);

  //! 重置所有累积数据
  void reset();

  //! 获取指定网格单元的通量
  //! \param cell_index 网格单元索引
  //! \return 通量值(如果该单元没有统计则返回0)
  double get_flux(int cell_index) const;

  //! 获取网格指针
  std::shared_ptr<SharedMeshGrid> get_grid() const { return grid_; }

private:
  //! 将空间位置转换为网格单元索引
  //! \param position 粒子位置 [x, y, z]
  //! \return 网格单元索引,如果超出范围返回-1
  int position_to_index(const std::array<double, 3>& position) const;

  //! 将一维索引转换为三维网格索引
  //! \param index 一维索引
  //! \return 三维网格索引 [ix, iy, iz]
  std::array<int, 3> index_to_grid(int index) const;

  //! 将三维网格索引转换为一维索引
  //! \param ix x方向索引
  //! \param iy y方向索引
  //! \param iz z方向索引
  //! \return 一维索引
  int grid_to_index(int ix, int iy, int iz) const;

  //! 输出到 HDF5 文件
  //! \param filename HDF5 文件名
  //! \param n_batches 总批次数
  void write_hdf5(const std::string& filename, int n_batches) const;

  //! 合并所有线程的局部数据到批次累积
  void merge_thread_local_data();

  // 数据成员
  std::shared_ptr<SharedMeshGrid> grid_; //!< 共享的网格对象

  // 当前批次的累积数据(稀疏存储)
  std::vector<double> make_zero_group_vector() const;
  int determine_group(double energy_eV, int mg_group) const;

  // 能群配置
  std::vector<double> energy_edges_;
  int n_groups_ {1};

  // 当前批次的累积数据(稀疏存储)
  std::unordered_map<int, std::vector<double>>
    batch_flux_group_; //!< 当前批次分群通量

  // 线程局部缓冲区(避免多线程竞争)
  std::vector<std::unordered_map<int, std::vector<double>>>
    thread_flux_group_; //!< 每个线程的局部分群累积

  // 跨批次聚合（分群）
  std::unordered_map<int, std::vector<double>>
    flux_group_sum_; //!< 各群通量总和
  std::unordered_map<int, std::vector<double>> flux_group_sum_sq_;

  // 兼容单群输出的标量统计
  std::unordered_map<int, double> flux_sum_;    //!< 通量总和 Σφ
  std::unordered_map<int, double> flux_sum_sq_; //!< 通量平方和 Σφ²

  int n_accumulated_batches_; //!< 已累积的批次数
};

} // namespace openmc

#endif // OPENMC_FLUX_MESH_H
