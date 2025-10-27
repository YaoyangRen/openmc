# 裂变矩阵运行控制指南

## 概述

裂变矩阵功能可以在特征值计算的不同阶段运行。本指南说明如何控制裂变矩阵在**活跃代**或**非活跃代**的运行。

## 当前实现状态

裂变矩阵功能已集成到 OpenMC 仿真流程中，目前通过 **默认仅在活跃代运行**。

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
- `n_realizations`: 统计的批次数
- `total_fissions`: 总裂变事件数
- `total_sources`: 总源粒子数

---

## 性能考虑

### 内存使用

裂变矩阵尺寸 = `n_cells × n_cells`

例如，网格 `100×100×100`：
- `n_cells = 1,000,000`
- 矩阵元素数 = `10^12`
- 内存需求 ≈ **8 TB** (double 类型)

⚠️ **建议**: 使用粗网格（1-5 cm）避免内存溢出

### 计算开销

数据记录操作：
- **源出生记录**: O(1) - 哈希表插入
- **裂变事件记录**: O(1) - 直接数组访问（有互斥锁）
- **批次管理**: O(n_cells²) - 矩阵累加

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

运行仿真时，控制台会输出：

### 初始化信息
```
FissionMatrix initialized:
  Bounds: [-10,-10,-10] to [10,10,10]
  Pitch: 1 cm
  Shape: [21,21,21]
  Total cells: 9261
  Matrix size: 9261 x 9261 = 85766121 elements
  Estimated memory: 651.32 MB
```

### 最终统计
```
Fission matrix written to: fission_matrix.h5
Statistics:
  Total fissions recorded: 1234567
  Total sources recorded: 50000
  Matrix size: 9261 x 9261
  Non-zero elements: 45123
  Sparsity: 99.47%
  Max normalized element: 0.0234
  Sum of normalized matrix: 1.982
```

---

## 常见问题

### Q1: 如何知道当前使用的是哪种模式？

查看 `src/simulation.cpp` 第 428 行附近，未被注释的代码即为当前模式。

### Q2: 可以在运行时切换模式吗？

不可以。模式在编译时确定，需要重新编译才能切换。

### Q3: 为什么推荐仅在活跃代运行？

活跃代的源分布已收敛，统计结果更准确。非活跃代的源分布仍在演化，可能引入偏差。

### Q4: 如何减少内存使用？

1. 增大网格间距（`pitch`）
2. 使用手动边界缩小计算区域
3. 仅在活跃代运行

---

## 相关文件

- **实现**: `src/fission_matrix.cpp`
- **头文件**: `include/openmc/fission_matrix.h`
- **集成**: `src/simulation.cpp`, `src/particle.cpp`, `src/physics.cpp`
- **文档**: `FISSION_MATRIX_CONTROL_GUIDE.md` (本文件)

---

## 更新日期

2025年10月27日
