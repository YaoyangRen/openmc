# 通量网格系统使用指南

# Flux Mesh System User Guide

## 概述 / Overview

通量网格系统用于统计空间通量分布,使用与裂变矩阵和传递函数完全相同的统一网格。

The flux mesh system tallies spatial flux distribution using the same unified grid as fission matrix and transfer function.

## 特性 / Features

1. **统一网格** - 与裂变矩阵和传递函数使用相同的网格,确保空间对齐
2. **稀疏存储** - 只存储非零通量单元,节省内存
3. **统计跟踪** - 自动计算均值和标准差
4. **径迹长度估计器** - 使用 Φ = w × d 进行通量估计
5. **HDF5输出** - 同时提供稀疏和密集格式

## 启用方法 / How to Enable

### 方法 1: C++ API

```cpp
// 在初始化代码中设置
settings::flux_mesh_on = true;
```

### 方法 2: Python API (openmc)

```python
import openmc

# 创建设置对象
settings = openmc.Settings()

# 添加自定义设置 (需要修改 Python 接口以支持)
# 或直接在 C++ 中启用
```

## 输出文件 / Output Files

### flux_mesh.h5

HDF5 文件包含以下数据集:

**网格参数 / Grid Parameters:**

- `grid_shape`: 网格尺寸 [nx, ny, nz]
- `grid_lower_left`: 下边界 [x, y, z] (cm)
- `grid_upper_right`: 上边界 [x, y, z] (cm)
- `grid_pitch`: 网格间距 [dx, dy, dz] (cm)
- `n_cells`: 总单元数
- `n_batches`: 批次数

**稀疏格式 / Sparse Format:**

- `cell_indices`: 非零单元索引列表
- `flux_mean`: 对应的平均通量 (n/cm²/source)
- `flux_std`: 对应的标准差 (n/cm²/source)

**密集格式 / Dense Format:**

- `flux_mean_dense`: 所有单元的平均通量 (0表示无统计)
- `flux_std_dense`: 所有单元的标准差

## 可视化 / Visualization

使用提供的 Python 脚本可视化通量分布:

```bash
python visualize_flux_mesh.py
```

### 生成的图像 / Generated Plots

1. **XY 切片** - `flux_xy_slice_z{index}.png`
   - 显示 Z 方向某一层的通量分布

2. **XZ 切片** - `flux_xz_slice_y{index}.png`
   - 显示 Y 方向某一层的通量分布

3. **相对误差** - `flux_error_z{index}.png`
   - 显示相对误差分布 (%)

4. **3D 散点图** - `flux_3d_scatter.png` (可选)
   - 三维可视化非零通量单元

## 代码集成 / Code Integration

### 文件清单 / File List

**新增文件 / New Files:**

- `include/openmc/flux_mesh.h` - FluxMesh 类定义
- `src/flux_mesh.cpp` - FluxMesh 实现
- `visualize_flux_mesh.py` - 可视化脚本

**修改文件 / Modified Files:**

- `include/openmc/simulation.h` - 添加 flux_mesh 声明
- `src/simulation.cpp` - 初始化和最终化 flux_mesh
- `include/openmc/settings.h` - 添加 flux_mesh_on 设置
- `src/settings.cpp` - 定义 flux_mesh_on 变量
- `src/tallies/tally_scoring.cpp` - 在径迹评分中累积通量
- `CMakeLists.txt` - 添加 flux_mesh.cpp 到编译

### 关键代码位置 / Key Code Locations

**通量累积 / Flux Accumulation:**

```cpp
// src/tallies/tally_scoring.cpp::score_tracklength_tally()
if (settings::flux_mesh_on && simulation::flux_mesh) {
  std::array<double, 3> position = {p.r().x, p.r().y, p.r().z};
  simulation::flux_mesh->accumulate(position, p.wgt(), distance);
}
```

**初始化 / Initialization:**

```cpp
// src/simulation.cpp::initialize_batch()
if (!simulation::flux_mesh && settings::flux_mesh_on) {
  simulation::flux_mesh = std::make_unique<FluxMesh>(shared_grid);
}
```

**批次处理 / Batch Processing:**

```cpp
// src/simulation.cpp::finalize_batch()
if (simulation::flux_mesh && settings::flux_mesh_on) {
  simulation::flux_mesh->end_batch(simulation::current_batch);
}
```

**最终化 / Finalization:**

```cpp
// src/simulation.cpp::finalize_batch() - 最后一批次
if (simulation::flux_mesh && settings::flux_mesh_on) {
  simulation::flux_mesh->finalize(settings::n_batches);
}
```

## 与其他系统的关系 / Relationship with Other Systems

```
SharedMeshGrid (统一网格)
    ├── FissionMatrix (裂变矩阵)
    ├── GreenFunctionMesh (传递函数)
    └── FluxMesh (通量分布) ← 新增
```

所有三个系统共享相同的网格参数,确保:

- 空间网格完全对齐
- 伴随源计算一致性
- 结果可以直接比较

## 性能考虑 / Performance Considerations

1. **内存使用** - 稀疏存储仅保存非零单元,内存占用与实际通量分布区域成正比
2. **计算开销** - 每个径迹都会累积通量,建议在需要时启用
3. **MPI 支持** - 目前简化实现,完整 MPI 归约需要进一步开发

## 示例工作流 / Example Workflow

```bash
# 1. 启用通量网格 (修改 settings.cpp 或在代码中设置)
settings::flux_mesh_on = true;

# 2. 编译项目
cd build
mingw32-make -j24

# 3. 运行模拟
./bin/openmc

# 4. 查看输出
# - flux_mesh.h5: 通量分布数据

# 5. 可视化结果
python ../visualize_flux_mesh.py

# 6. 查看生成的图像
# - flux_xy_slice_z*.png
# - flux_xz_slice_y*.png
# - flux_error_z*.png
```

## 故障排除 / Troubleshooting

**问题 1: 编译错误**

- 确保 flux_mesh.cpp 已添加到 CMakeLists.txt
- 检查所有头文件路径正确

**问题 2: 运行时未生成文件**

- 检查 `settings::flux_mesh_on` 是否设置为 `true`
- 查看控制台输出是否有 "通量网格已初始化" 消息

**问题 3: HDF5 文件为空**

- 确保模拟运行到最后一批次
- 检查是否有粒子经过网格区域

**问题 4: 可视化脚本报错**

- 安装必需的 Python 包: `pip install h5py numpy matplotlib`
- 确认 flux_mesh.h5 文件存在

## 未来改进 / Future Improvements

1. **完整 MPI 支持** - 实现跨进程的稀疏数据归约
2. **能量分组** - 支持多群通量统计
3. **在线可视化** - 批次间实时查看通量分布
4. **自适应网格** - 根据通量梯度动态调整网格密度
5. **Python API** - 添加 Python 接口控制开关

---

**版本**: 1.0  
**日期**: 2024  
**作者**: OpenMC Flux Mesh Development Team
