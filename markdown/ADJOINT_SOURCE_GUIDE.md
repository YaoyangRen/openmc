# 伴随源迭代计算功能

**版本**: v1.0  
**日期**: 2025年11月4日

---

## 📖 理论基础

### 伴随源分布方程

伴随裂变源分布 $\vec{I^*}$ 满足特征值方程：

$$\vec{I^*} = \frac{1}{k} \vec{F^*} \vec{I^*} = \frac{1}{k} \vec{F^T} \vec{I^*}$$

其中：

- $\vec{I^*}$ 是伴随源分布向量（重要性函数）
- $\vec{F^T}$ 是裂变矩阵的转置
- $k$ 是有效增殖因子（主特征值）
- $\vec{F^*}$ 表示伴随裂变矩阵

### 物理意义

- **正向问题**：$\vec{S} = \frac{1}{k} \vec{F} \vec{S}$  
  计算中子在空间的分布

- **伴随问题**：$\vec{I^*} = \frac{1}{k} \vec{F^T} \vec{I^*}$  
  计算每个空间位置对系统的"重要性"

伴随源 $I^*_i$ 表示在位置 $i$ 产生一个中子对整个系统的贡献。

---

## 🚀 功能特性

### 1. 幂迭代法求解

采用经典幂迭代（Power Iteration）算法：

```
初始化: I* = I_0
重复直到收敛:
    1. I_new = F^T × I*
    2. k = sum(I_new)
    3. I* = I_new / k
    4. 检查 |k - k_old| < tolerance
```

### 2. 初始值设置

支持两种初始化方式：

#### (1) 均匀分布 (`uniform`)

```
I*_i = 1,  ∀i
```

- **优点**: 稳健、无偏、适合所有问题
- **缺点**: 可能需要更多迭代次数
- **推荐**: 默认选项

#### (2) 正向源分布 (`forward`)

```
I*_i = S_i / Σ S_i
```

- **优点**: 收敛速度快（如果正向源已收敛）
- **缺点**: 依赖正向源质量
- **推荐**: 正向源已收敛时使用

### 3. Batch级伴随源迭代（新特性）

**在每个batch结束后立即更新伴随源分布**

#### 启用方式

```cpp
// 启用batch级迭代：每个batch执行10次迭代
fission_matrix->enable_batch_adjoint_iteration(
    true,   // 启用
    10,     // 每batch迭代次数
    1.0e-6  // 收敛容差
);
```

#### 工作原理

- 每完成一个batch，立即使用累积的裂变矩阵进行伴随源迭代
- 伴随源随着统计的积累逐步收敛
- 记录每个batch后的 $k_{adjoint}$ 值，可分析收敛过程

#### 优势

- **实时更新**: 伴随源与裂变矩阵同步演化
- **收敛监控**: 可观察伴随源如何随batch收敛
- **灵活控制**: 可调整每batch迭代次数平衡精度和性能

### 4. 自适应收敛判据

- **默认容差**: $10^{-6}$
- **最大迭代次数**: 1000
- **收敛判据**: $|k_{\text{new}} - k_{\text{old}}| < \text{tolerance}$

---

## 💻 使用方法

### 方法1: Batch级迭代（推荐）

```cpp
#include "openmc/fission_matrix.h"

// 1. 创建裂变矩阵对象
auto fission_matrix = std::make_unique<FissionMatrix>(
    0.5,        // pitch: 网格间距 (cm)
    100,        // max_batches: 最大批次数
    true        // auto_bounds: 自动边界检测
);

// 2. 启用batch级伴随源迭代
fission_matrix->enable_batch_adjoint_iteration(
    true,   // 启用
    10,     // 每batch迭代10次
    1.0e-6  // 收敛容差
);

// 3. 在模拟过程中记录裂变事件
// (自动调用 record_source_birth 和 record_fission_event)
// 每个batch结束时自动执行伴随源迭代

// 4. 最终化并保存（可选：额外的精细化迭代）
fission_matrix->compute_adjoint_source("uniform", 100, 1e-8);
fission_matrix->finalize("fission_matrix.h5");
```

### 方法2: 传统方法（仅最终化时计算）

```cpp
#include "openmc/fission_matrix.h"

// 1. 创建裂变矩阵对象
auto fission_matrix = std::make_unique<FissionMatrix>(
    0.5,        // pitch: 网格间距 (cm)
    100,        // max_batches: 最大批次数
    true        // auto_bounds: 自动边界检测
);

// 2. 在模拟过程中记录裂变事件
// (自动调用 record_source_birth 和 record_fission_event)

// 3. 计算伴随源分布（仅在最后执行）
fission_matrix->compute_adjoint_source(
    "uniform",  // 初始猜测: "uniform" 或 "forward"
    1000,       // 最大迭代次数
    1.0e-6      // 收敛容差
);

// 4. 最终化并保存
fission_matrix->finalize("fission_matrix.h5");
```

### 完整示例

```cpp
// simulation.cpp 中的调用示例

// 初始化阶段
if (settings::run_fission_matrix) {
    simulation::fission_matrix = std::make_unique<FissionMatrix>(
        settings::fission_matrix_pitch,
        settings::n_batches
    );
}

// 批次循环中
for (int batch = 0; batch < settings::n_batches; ++batch) {
    simulation::fission_matrix->start_new_batch(batch);
    
    // ... 粒子输运 ...
    
    // 在每个源粒子初始化时
    simulation::fission_matrix->record_source_birth(
        particle.r(), particle.id()
    );
    
    // 在裂变事件时
    simulation::fission_matrix->record_fission_event(
        particle.r(), nu_fission, particle.source_particle_id()
    );
}

// 模拟结束后
// 计算伴随源（使用均匀初始化）
simulation::fission_matrix->compute_adjoint_source("uniform", 1000, 1e-6);

// 或使用正向源初始化（更快收敛）
// simulation::fission_matrix->compute_adjoint_source("forward", 1000, 1e-6);

simulation::fission_matrix->finalize("fission_matrix.h5");
```

---

## 📊 输出数据

### HDF5 文件结构

```text
fission_matrix.h5
├── row_indices          [int array]    - COO格式行索引
├── col_indices          [int array]    - COO格式列索引
├── data_raw             [double array] - 原始裂变矩阵数据
├── data_normalized      [double array] - 归一化裂变矩阵
├── source_counts        [double array] - 每个单元的源粒子数
├── adjoint_source       [double array] - 伴随源分布 I*
├── k_adjoint_history    [double array] - batch级k_adjoint历史 (可选)
│
└── 属性 (Attributes):
    ├── filetype         = "fission_matrix_sparse"
    ├── version          = "2.0"
    ├── storage_format   = "COO"
    ├── n_cells          - 总单元数
    ├── nnz              - 非零元素数
    ├── k_adjoint        - 伴随k值
    ├── adjoint_iterations - 总迭代次数
    ├── adjoint_converged  - 是否收敛
    ├── batch_adjoint_enabled   - 是否启用batch级迭代 (可选)
    └── adjoint_iter_per_batch  - 每batch迭代次数 (可选)
```

### Python 读取示例

```python
import h5py
import numpy as np
from scipy.sparse import coo_matrix

# 读取伴随源
with h5py.File('fission_matrix.h5', 'r') as f:
    # 读取伴随源分布
    adjoint_source = f['adjoint_source'][:]
    k_adjoint = f.attrs['k_adjoint']
    iterations = f.attrs['adjoint_iterations']
    
    # 读取网格信息
    shape = f['shape'][:]
    
    # 重塑为3D数组
    I_star_3d = adjoint_source.reshape(shape)
    
    print(f"k_adjoint = {k_adjoint:.8f}")
    print(f"Iterations = {iterations}")
    print(f"Max importance = {adjoint_source.max():.6e}")
```

---

## 📈 可视化分析

### 使用提供的Python脚本

```bash
python analyze_adjoint_source.py fission_matrix.h5
```

**功能**:

1. 读取并显示伴随源统计信息
2. 生成2D切片热图
3. 比较正向源和伴随源分布
4. **绘制batch级迭代收敛历史** (新增)
5. 分析收敛历史
6. 计算不同初始化方式的收敛性能

**输出**:

- `adjoint_source_distribution.png` - 伴随源分布图
- `forward_adjoint_comparison.png` - 正向/伴随对比图
- `batch_adjoint_convergence.png` - batch级收敛历史图 (如果启用了batch级迭代)

---

## 🔬 算法详解

### 伴随矩阵-向量乘法

对于稀疏矩阵 $F$ (COO格式)：

```
F[i][j] 存储为 key = i × n_cells + j

计算 F^T × I*:
    for each (key, F_ij) in sparse_matrix:
        i = key / n_cells  # 源单元
        j = key % n_cells  # 裂变单元
        
        # F^T[j][i] = F[i][j]
        result[j] += F_ij × I*[i]
```

**时间复杂度**: $O(\text{nnz})$，其中 nnz 是非零元素数

### 收敛性分析

**理论保证**:

- 幂迭代法对主特征值收敛
- 收敛率: $\lambda_1 / \lambda_2$（特征值比）
- 典型问题: 50-200 次迭代收敛

**加速技巧**:

1. 使用正向源初始化（减少50%迭代）
2. 提高容差到 $10^{-8}$（更精确）
3. 预条件（未实现）

---

## 📋 应用场景

### 1. 扰动分析

计算参数扰动对 $k_{eff}$ 的影响：

$$\delta k = \frac{\langle I^*, \delta F S \rangle}{\langle I^*, F S \rangle}$$

### 2. 重要性采样

- 按 $I^*$ 分布采样源位置
- 提高关键区域的统计精度
- 减少方差

### 3. 探测器响应优化

- 计算探测器位置的最佳布置
- $I^*$ 高的区域对探测器影响大

### 4. 燃料装载优化

- 识别对 $k_{eff}$ 影响最大的区域
- 优化燃料布置方案

---

## ⚠️ 注意事项

### 1. 裂变矩阵质量

- **建议**: 在数个非活跃代后开始统计裂变矩阵
- **原因**: 避免初始源分布偏差
- **典型设置**: 跳过前 20-50 个非活跃代

### 2. 内存使用

- **稀疏存储**: 仅存储非零元素
- **典型问题**: 99% 稀疏度
- **示例**: 50³ 网格仅需 ~36 MB（而非 117 GB）

### 3. 收敛判据

```cpp
// 严格收敛
compute_adjoint_source("uniform", 2000, 1e-8);

// 快速估算
compute_adjoint_source("forward", 500, 1e-5);
```

### 4. 调试输出

每50次迭代输出进度：

```
Iteration   50: k_adj = 1.023456, dk = 3.45e-04
Iteration  100: k_adj = 1.023489, dk = 8.76e-05
...
Converged at iteration 287
  Final k_adjoint = 1.02349876
  Final dk = 9.87e-07
```

---

## 🧪 验证测试

### 单元测试

```cpp
// 验证伴随源归一化
double sum = 0.0;
for (auto val : adjoint_source_) {
    sum += val;
}
assert(std::abs(sum - 1.0) < 1e-10);

// 验证k值一致性（正向与伴随）
assert(std::abs(k_adjoint - k_forward) < 0.01);
```

### 与解析解对比

对于简单几何（无限均匀介质）:

- 伴随源应为常数
- $I^*_i = 1/N$ for all $i$

---

## 🔧 故障排除

### 问题 1: 未收敛

**症状**: "Maximum iterations reached without convergence"

**解决**:

```cpp
// 增加迭代次数
compute_adjoint_source("uniform", 5000, 1e-6);

// 或放宽容差
compute_adjoint_source("uniform", 1000, 1e-5);
```

### 问题 2: k值异常

**症状**: k_adjoint 与 k_forward 差异很大

**原因**: 裂变矩阵统计不充分

**解决**: 增加批次数或每批次源粒子数

### 问题 3: 伴随源为零

**症状**: adjoint_source 全为零

**原因**: 裂变矩阵为空

**检查**:

```cpp
if (fission_matrix_sparse_.empty()) {
    std::cerr << "Error: No fission events recorded!" << std::endl;
}
```

---

## 📞 获取帮助

如有问题，请检查：

1. 裂变矩阵是否正确统计（nnz > 0）
2. 初始猜测是否合理
3. 收敛容差是否过严
4. 迭代次数是否足够

**调试信息输出**: 设置 `verbosity >= 7` 查看详细日志

---

**文档版本**: v1.0  
**最后更新**: 2025年11月4日  
**维护者**: RYY
