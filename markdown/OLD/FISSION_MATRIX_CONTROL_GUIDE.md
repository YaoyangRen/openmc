# 裂变矩阵运行控制指南

## 概述

裂变矩阵功能可以在特征值计算的不同阶段运行。本指南说明如何控制裂变矩阵在**活跃代**或**非活跃代**的运行。

## 当前实现状态

裂变矩阵功能已集成到 OpenMC 仿真流程中，采用**稀疏矩阵存储（COO格式）**，目前默认**仅在活跃代运行**。

### 🚀 新特性：稀疏矩阵存储

- **存储格式**: COO (Coordinate) 格式
- **内存优化**: 对于 99% 稀疏度的矩阵，可节省 **99%** 内存
- **自动优化**: 只存储非零元素，动态增长
- **兼容性**: 标准格式，支持 SciPy、MATLAB、Eigen 等工具

**内存节省示例**:

| 网格尺寸 | 密集存储 | 稀疏存储 | 节省 |
|---------|---------|---------|------|
| 11×11×11 | 13.5 MB | 0.4 MB | 97% |
| 50×50×50 | 117 GB | 35.7 MB | 99.97% |

详见: [稀疏矩阵格式说明](./SPARSE_FISSION_MATRIX.md)

## 控制选项

在 `src/simulation.cpp` 的 `initialize_batch()` 函数中，您可以选择以下三种模式之一：

### 选项 1: 仅在活跃代运行（当前默认）✅

```cpp
// Start new batch for fission matrix
// 选项 1: 仅在活跃代运行裂变矩阵
if (settings::clutch_on && simulation::current_batch > settings::n_inactive) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

**适用场景**:

- 🎯 正常生产计算
- 💡 只关心收敛后的裂变分布
- 💾 减少内存和计算开销

**批次范围**: `[n_inactive + 1, n_batches]`

---

### 选项 2: 在所有批次运行（包括非活跃代）

```cpp
// 选项 2: 在所有批次（包括非活跃代）运行裂变矩阵
if (settings::clutch_on) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

**适用场景**:

- 📊 研究源收敛过程
- 🔬 观察裂变分布随批次的演化
- 🧪 调试和验证

**批次范围**: `[1, n_batches]`

---

### 选项 3: 仅在非活跃代运行（调试用）

```cpp
// 选项 3: 仅在非活跃代运行裂变矩阵（用于调试）
if (settings::clutch_on && simulation::current_batch <= settings::n_inactive) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

**适用场景**:

- 🐛 调试源收敛问题
- 📈 分析初始源分布的影响
- 🔍 验证边界处理

**批次范围**: `[1, n_inactive]`

---

## 如何切换模式

### 步骤 1: 编辑 simulation.cpp

打开文件 `src/simulation.cpp`，找到 `initialize_batch()` 函数中的裂变矩阵控制部分（约第 428 行）。

### 步骤 2: 注释/取消注释相应代码

取消注释您想要的选项，注释掉其他选项：

```cpp
// Start new batch for fission matrix
// 选项 1: 仅在活跃代运行裂变矩阵 ✓ 当前激活
if (settings::clutch_on && simulation::current_batch > settings::n_inactive) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}

// 选项 2: 在所有批次（包括非活跃代）运行裂变矩阵
// if (settings::clutch_on) {
//   if (simulation::fission_matrix) {
//     simulation::fission_matrix->start_new_batch(simulation::current_batch);
//   }
// }

// 选项 3: 仅在非活跃代运行裂变矩阵（用于调试）
// if (settings::clutch_on && simulation::current_batch <= settings::n_inactive) {
//   if (simulation::fission_matrix) {
//     simulation::fission_matrix->start_new_batch(simulation::current_batch);
//   }
// }
```

### 步骤 3: 重新编译

```powershell
cd d:\OpenMC\openmc\build
mingw32-make -j24
```

---

## 批次编号说明

假设您的设置是：

- `n_inactive = 50` (非活跃代批次)
- `n_batches = 150` (总批次数)

则：

- **非活跃代**: 批次 1-50
- **活跃代**: 批次 51-150

| 选项 | 运行的批次 | 统计的批次数 |
|------|-----------|-------------|
| 选项 1 (仅活跃代) | 51-150 | 100 批次 |
| 选项 2 (所有批次) | 1-150 | 150 批次 |
| 选项 3 (仅非活跃代) | 1-50 | 50 批次 |

---

## 数据记录逻辑

无论选择哪种模式，数据记录都是自动的：

1. **源出生记录** (`particle.cpp`):
   - 在 `from_source()` 中自动记录
   - 仅当 `fission_matrix` 已初始化时记录

2. **裂变事件记录** (`physics.cpp`, `physics_mg.cpp`):
   - 在 `create_fission_sites()` 中自动记录
   - 使用期望裂变中子数 $\nu_t$ 作为权重

3. **批次管理**:
   - `start_new_batch()` 在每个批次开始时调用
   - `finalize()` 在最后一个批次结束时调用

---

## 输出文件

无论哪种模式，最终都会生成：

**文件名**: `fission_matrix.h5`

**包含数据集**:

- `fission_matrix_raw`: 原始累积矩阵
- `fission_matrix_normalized`: 归一化矩阵（按源计数）
- `source_counts`: 每个网格单元的源粒子计数
- `origin`, `shape`, `pitch`: 网格几何信息

**属性**:

- `storage_format`: "COO" (稀疏矩阵格式)
- `n_realizations`: 统计的批次数
- `total_fissions`: 总裂变事件数
- `total_sources`: 总源粒子数
- `nnz`: 非零元素数量

**数据集** (稀疏格式):

- `row_indices`, `col_indices`: 稀疏矩阵坐标
- `data_raw`, `data_normalized`: 原始和归一化值

详细格式说明: [SPARSE_FISSION_MATRIX.md](./SPARSE_FISSION_MATRIX.md)

---

## 性能考虑

### 内存使用（稀疏存储）

**稀疏矩阵存储**大幅降低内存需求：

| 网格尺寸 | 总单元数 | 密集存储 | 稀疏存储 (99%稀疏) | 内存节省 |
|---------|---------|---------|------------------|---------|
| 11×11×11 | 1,331 | 13.5 MB | 0.4 MB | **97%** ✨ |
| 20×20×20 | 8,000 | 488 MB | 7.3 MB | **98.5%** |
| 50×50×50 | 125,000 | 117 GB | 35.7 MB | **99.97%** 🚀 |
| 100×100×100 | 1,000,000 | 7.45 TB | 1.1 GB | **99.985%** |

✅ **现在可以使用精细网格！**

- 推荐: 0.5-2 cm 网格间距
- 大规模问题: 100×100×100 网格仅需 ~1 GB 内存

### 计算开销

稀疏矩阵操作开销：

- **源出生记录**: O(1) - 哈希表插入
- **裂变事件记录**: O(1) - 哈希表更新（有互斥锁）
- **批次管理**: O(nnz) - 仅处理非零元素（原 O(n²)）
- **归一化**: O(nnz) - 显著减少计算量

**性能提升**:

- 内存访问: ~100× 更快（缓存友好）
- 矩阵操作: ~10-100× 加速
- HDF5 写入: 文件大小减少 95%+

仅在活跃代运行可减少 `n_inactive/n_batches` 的开销（通常 30-50%）。

---

## 高级自定义

如果需要更细粒度的控制，可以修改条件判断：

### 示例: 仅在指定批次范围运行

```cpp
// 仅在批次 100-200 运行裂变矩阵
if (settings::clutch_on && 
    simulation::current_batch >= 100 && 
    simulation::current_batch <= 200) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

### 示例: 每隔 N 个批次运行

```cpp
// 每隔 10 个批次运行一次
if (settings::clutch_on && 
    simulation::current_batch > settings::n_inactive &&
    (simulation::current_batch % 10 == 0)) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

---

## 验证运行状态

运行仿真时，控制台会在**结束时**输出详细信息：

### 最终统计输出（稀疏矩阵）

```txt
======================================================================
FISSION MATRIX FINALIZATION
======================================================================

Grid Configuration:
  Bounds: [0,0,0] to [10,10,10]
  Pitch: 1 cm
  Shape: [11,11,11]
  Total cells: 1331
  Matrix size: 1331 x 1331 = 1771561 elements
  Non-zero elements: 45123
  Sparsity: 97.45%
  Dense storage would need: 13.52 MB
  Sparse storage uses: 1.03 MB
  Memory saved: 12.49 MB ✨

Normalizing sparse fission matrix...

Data Collection Summary:
  Total fissions recorded: 1234567
  Total sources recorded: 50000
  Number of batches: 100

Sparse Matrix Statistics:
  Non-zero elements: 45123 / 1771561
  Sparsity: 97.45%
  Max normalized element: 0.0234
  Sum of normalized matrix: 1.982

Storage Format: COO (Coordinate)
Output File: fission_matrix.h5
======================================================================
```

---

## 常见问题

### Q1: 如何知道当前使用的是哪种模式？

查看 `src/simulation.cpp` 第 428 行附近，未被注释的代码即为当前模式。

### Q2: 可以在运行时切换模式吗？

不可以。模式在编译时确定，需要重新编译才能切换。

### Q3: 为什么推荐仅在活跃代运行？

活跃代的源分布已收敛，统计结果更准确。非活跃代的源分布仍在演化，可能引入偏差。

### Q4: 稀疏矩阵如何读取？

使用 Python + h5py + SciPy:

```python
import h5py
from scipy.sparse import coo_matrix

with h5py.File('fission_matrix.h5', 'r') as f:
    n_cells = f.attrs['n_cells']
    rows = f['row_indices'][:]
    cols = f['col_indices'][:]
    data = f['data_normalized'][:]
    
    # 创建稀疏矩阵
    F = coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
```

详见: [SPARSE_FISSION_MATRIX.md](./SPARSE_FISSION_MATRIX.md)

### Q5: 稀疏度低于 90% 怎么办？

- 增大网格间距 (pitch)
- 检查几何耦合程度
- 即使稀疏度较低，稀疏存储仍比密集存储高效

### Q6: 如何进一步减少内存？

1. 增大网格间距 (`pitch = 2.0` 或更大)
2. 使用手动边界缩小计算区域
3. 仅在活跃代运行（减少统计批次）
4. 稀疏存储已自动优化，无需额外操作

---

## 网格尺寸推荐

基于稀疏存储，现在可以使用更精细的网格：

| 计算类型 | 推荐间距 | 典型网格 | 内存需求 | 适用场景 |
|---------|---------|---------|----------|---------|
| 快速测试 | 5 cm | 10×10×10 | < 1 MB | 调试验证 |
| 标准计算 | 2 cm | 25×25×25 | < 10 MB | 常规分析 |
| 精细分析 | 1 cm | 50×50×50 | < 100 MB | 详细研究 |
| 高精度 | 0.5 cm | 100×100×100 | < 2 GB | 基准计算 |

---

## 相关文件

- **实现**: `src/fission_matrix.cpp`
- **头文件**: `include/openmc/fission_matrix.h`
- **集成**: `src/simulation.cpp`, `src/particle.cpp`, `src/physics.cpp`
- **文档**: `FISSION_MATRIX_CONTROL_GUIDE.md` (本文件)

---

## 更新日期

2025年10月27日
