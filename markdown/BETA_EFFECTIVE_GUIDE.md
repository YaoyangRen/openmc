# 有效缓发中子份额 (β_eff) 计算模块使用手册

## 目录

- [概述](#概述)
- [理论背景](#理论背景)
- [使用方法](#使用方法)
- [输出文件格式](#输出文件格式)
- [开发路线图](#开发路线图)
- [示例与验证](#示例与验证)
- [故障排除](#故障排除)

---

## 概述

### 功能简介

β_eff 计算模块用于计算堆芯的**有效缓发中子份额**，这是反应堆动力学中的关键参数，直接影响反应堆的控制性能和安全特性。

### 物理意义

- **β_i,eff**: 第 i 组缓发中子的有效份额，考虑了空间和能量的共轭通量加权
- **β_total**: 总有效缓发中子份额，用于点堆动力学方程
- **重要性**: β_eff 越大，反应堆越容易控制；β_eff 越小，反应性引入的风险越高

### 当前状态

**Phase 1 (已实现)**: 基于固定核数据的单能群计算

- 核数据: U-235 热中子裂变参数
- 能群: 单能群近似 (空间积分)
- 不确定度: 暂无统计传播

---

## 理论背景

### 基本公式

有效缓发中子份额的定义基于微扰理论:

$$
\beta_{i,\text{eff}} = \frac{\int_V \phi^*(r) \cdot \chi_{d,i} \cdot \nu_{d,i} \cdot \Sigma_f(r) \cdot \phi(r) \, dV}{\int_V \phi^*(r) \cdot \chi_p \cdot \nu_{\text{total}} \cdot \Sigma_f(r) \cdot \phi(r) \, dV}
$$

总有效缓发中子份额:

$$
\beta_{\text{eff}} = \sum_{i=1}^{6} \beta_{i,\text{eff}}
$$

### 符号说明

| 符号 | 含义 | Phase 1 取值 |
|------|------|-------------|
| $\phi(r)$ | 正向中子通量 | 从 `flux_mesh.h5` 读取 |
| $\phi^*(r)$ | 共轭中子通量 | 从 `adjoint_flux.h5` 读取 |
| $\nu_{\text{total}}$ | 总裂变中子数 | 2.43 (U-235) |
| $\nu_{d,i}$ | 第 i 组缓发中子数 | 见下表 |
| $\chi_p$ | 瞬发中子能谱 | 1.0 (归一化) |
| $\chi_{d,i}$ | 第 i 组缓发中子能谱 | 1.0 (归一化) |
| $\Sigma_f(r)$ | 宏观裂变截面 | 1.0 (相对值，约掉) |
| $V$ | 体积元 | $\text{pitch}^3$ |

### U-235 缓发中子数据

| 组别 | $\nu_{d,i}$ | 半衰期 (s) |
|------|-------------|-----------|
| 1 | 0.000215 | 55.72 |
| 2 | 0.001424 | 22.72 |
| 3 | 0.001274 | 6.22 |
| 4 | 0.002568 | 2.30 |
| 5 | 0.000748 | 0.614 |
| 6 | 0.000273 | 0.230 |
| **总和** | **0.00650** | - |

理论值: $\beta_{\text{theory}} = \sum \nu_{d,i} / \nu_{\text{total}} = 0.00650 / 2.43 \approx 0.00267$ (0.267%)

---

## 使用方法

### 前置条件

1. **开启 Clutch 模式** (在 `settings.xml` 中):

```xml
<clutch>
  <on>true</on>
</clutch>
```

2. **确保通量网格已配置** (自动开启):
   - `flux_mesh_on` 会随 `clutch_on` 自动启用
   - 网格参数: 19×19×19, pitch=1cm, 范围 [-8.791, 8.791]

3. **运行完整模拟**:
   - 至少运行到 `n_batches`，确保生成 `flux_mesh.h5` 和 `adjoint_flux.h5`

### 自动计算流程

β_eff 计算在模拟的最后一个 batch 自动触发，无需额外操作:

```
Batch 100/100 完成
  ↓
计算传递函数 (Green Function)
  ↓
计算共轭通量 (Adjoint Flux)
  ↓
【自动触发】计算 β_eff
  ↓
输出 beta_eff.h5
```

### 控制台输出示例

```
======================================================================
COMPUTING EFFECTIVE DELAYED NEUTRON FRACTION
======================================================================
Phase 1: Using fixed U-235 thermal neutron data
  Energy groups: 1 (single-group approximation)

Reading flux data from flux_mesh.h5...
  Grid: [19, 19, 19], pitch = 1 cm
  Non-zero cells: 2847

Reading adjoint flux data from adjoint_flux.h5...
  Grid: [19, 19, 19], pitch = 1 cm
  Non-zero cells: 2951

Grid validation: OK

Computing delayed neutron contributions...
  Group 1: β_1 = 8.8480e-05
  Group 2: β_2 = 5.8600e-04
  Group 3: β_3 = 5.2430e-04
  Group 4: β_4 = 1.0570e-03
  Group 5: β_5 = 3.0780e-04
  Group 6: β_6 = 1.1230e-04

β_eff computation completed successfully.
  β_eff = 0.00267

  Results written successfully
  Theoretical β = 0.00267
  Computed β_eff = 0.00267
  Relative difference = 0.12%
```

---

## 输出文件格式

### HDF5 结构

```
beta_eff.h5
├── beta_i [6]                    // 数据集: 各组 β_i,eff
├── beta_total [1]                // 数据集: 总 β_eff
├── uncertainty [6]               // 数据集: 不确定度 (Phase 1 全为 0)
│
├── metadata/                     // 组: 元数据
│   ├── flux_file                 // 属性: "flux_mesh.h5"
│   ├── adjoint_flux_file         // 属性: "adjoint_flux.h5"
│   ├── computation_time          // 属性: Unix 时间戳
│   ├── nuclear_data_source       // 属性: "U-235 thermal (Phase 1)"
│   ├── energy_groups             // 属性: 1
│   ├── grid_shape [3]            // 数据集: [19, 19, 19]
│   └── grid_pitch                // 属性: 1.0
│
└── diagnostics/                  // 组: 诊断信息
    ├── numerator [6]             // 数据集: 各组分子项
    ├── denominator [1]           // 数据集: 分母项
    ├── normalization_factor      // 属性: 体积元 (pitch³)
    ├── beta_theoretical          // 属性: 理论值 0.00267
    ├── relative_difference_percent // 属性: 相对偏差 %
    ├── nu_total                  // 属性: 2.43
    ├── nu_prompt                 // 属性: 2.42
    └── nu_delayed [6]            // 数据集: [0.000215, ...]
```

### Python 读取示例

```python
import h5py
import numpy as np

# 读取 beta_eff.h5
with h5py.File('beta_eff.h5', 'r') as f:
    # 主要结果
    beta_i = f['beta_i'][:]
    beta_total = f['beta_total'][0]
    uncertainty = f['uncertainty'][:]
    
    # 元数据
    flux_file = f['metadata'].attrs['flux_file']
    comp_time = f['metadata'].attrs['computation_time']
    grid_shape = f['metadata/grid_shape'][:]
    
    # 诊断信息
    numerators = f['diagnostics/numerator'][:]
    denominator = f['diagnostics/denominator'][0]
    beta_theory = f['diagnostics'].attrs['beta_theoretical']
    rel_diff = f['diagnostics'].attrs['relative_difference_percent']
    
print(f"β_eff = {beta_total:.5f}")
print(f"理论值 = {beta_theory:.5f}")
print(f"相对偏差 = {rel_diff:.2f}%")
print(f"\n各组贡献:")
for i, bi in enumerate(beta_i, 1):
    print(f"  β_{i} = {bi:.5e} ({bi/beta_total*100:.1f}%)")
```

### 可视化脚本

创建 `visualize_beta_eff.py`:

```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File('beta_eff.h5', 'r') as f:
    beta_i = f['beta_i'][:]
    beta_total = f['beta_total'][0]
    nu_delayed = f['diagnostics/nu_delayed'][:]

groups = np.arange(1, 7)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# 左图: β_i 分布
ax1.bar(groups, beta_i, color='steelblue', alpha=0.7)
ax1.axhline(beta_total/6, color='red', linestyle='--', 
            label=f'平均值 = {beta_total/6:.5e}')
ax1.set_xlabel('缓发中子组')
ax1.set_ylabel('β_i,eff')
ax1.set_title('各组有效缓发中子份额')
ax1.legend()
ax1.grid(axis='y', alpha=0.3)

# 右图: 与 ν_d 对比
ax2.scatter(nu_delayed, beta_i, s=100, alpha=0.6)
for i, (x, y) in enumerate(zip(nu_delayed, beta_i), 1):
    ax2.annotate(f'组{i}', (x, y), xytext=(5, 5), 
                 textcoords='offset points')
ax2.set_xlabel('ν_d,i (每次裂变缓发中子数)')
ax2.set_ylabel('β_i,eff')
ax2.set_title('β_i,eff vs ν_d,i 相关性')
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('beta_eff_analysis.png', dpi=300)
plt.show()

print(f"β_eff = {beta_total:.5f}")
print(f"β_total/Σν_d = {beta_total/nu_delayed.sum():.3f} (应接近 1.0)")
```

---

## 开发路线图

### Phase 1: 固定核数据单能群 ✅ (当前版本)

**特点:**

- U-235 热中子裂变固定参数
- 单能群空间积分
- χ_p = χ_d,i = 1.0 归一化假设
- 稀疏存储优化

**适用场景:**

- 快速原型验证
- U-235 主导的热中子堆
- 相对变化趋势分析

**限制:**

- 不适用于 Pu-239、MOX 燃料
- 忽略能谱效应
- 无统计不确定度

### Phase 2: 材料依赖核数据 🔄 (计划中)

**目标:**
从 OpenMC 材料库提取实际核数据

**实现要点:**

```cpp
// 伪代码示例
for (auto& cell : cells) {
  auto mat = cell.material;
  for (auto& nuclide : mat.nuclides) {
    if (nuclide.is_fissile()) {
      nu_total += nuclide.fraction * nuclide.nu_total();
      for (int i = 0; i < 6; i++) {
        nu_delayed[i] += nuclide.fraction * nuclide.nu_delayed(i);
      }
    }
  }
}
```

**优势:**

- 支持混合燃料 (UO₂, MOX, 钍燃料等)
- 考虑燃耗过程中的核素变化
- 更准确的物理模型

### Phase 3: 多能群扩展 🔄 (计划中)

**公式扩展:**

$$
\beta_{i,\text{eff}} = \frac{\sum_{g=1}^G \int_V \phi^*_g(r) \cdot \chi_{d,i,g'} \cdot \nu_{d,i} \cdot \Sigma_{f,g}(r) \cdot \phi_g(r) \, dV}{\sum_{g=1}^G \int_V \phi^*_g(r) \cdot \chi_{p,g'} \cdot \nu_{\text{total}} \cdot \Sigma_{f,g}(r) \cdot \phi_g(r) \, dV}
$$

**数据结构修改:**

```cpp
class BetaEffective {
  std::vector<std::vector<double>> beta_i_group_; // [6][n_groups]
  std::vector<double> beta_i_integrated_;         // [6]
  // ...
};
```

**输出增强:**

```
beta_eff.h5
├── beta_i_by_group [6, n_groups]  // 能群分辨结果
├── beta_i [6]                     // 能群积分结果
└── energy_grid [n_groups+1]       // 能量边界
```

### Phase 4: 统计不确定度传播 🔄 (计划中)

**方法:**

- Bootstrap 重采样
- 协方差传播
- 置信区间估计

**输出示例:**

```
β_eff = 0.00267 ± 0.00003 (1σ)
```

**HDF5 更新:**

```cpp
write_dataset(file_id, "uncertainty", uncertainty_vec); // 不再全为 0
write_dataset(diag_group, "covariance_matrix", cov_matrix);
```

---

## 示例与验证

### 典型输出值参考

| 堆型 | β_eff | 备注 |
|------|-------|------|
| 热中子堆 (U-235) | 0.0064-0.0068 | Phase 1 理论值 |
| 压水堆 (低富集度) | 0.0050-0.0060 | Phase 2+ 需要 |
| 快堆 (Pu-239) | 0.0020-0.0035 | Phase 2+ 需要 |
| 钍堆 (U-233) | 0.0026-0.0030 | Phase 2+ 需要 |

### 验证检查清单

**1. 数值合理性:**

```python
# 检查 β_eff 范围
assert 0.002 < beta_total < 0.010, "β_eff 超出物理范围"

# 检查各组非负
assert all(bi >= 0 for bi in beta_i), "存在负值"

# 检查和一致性
assert abs(sum(beta_i) - beta_total) < 1e-10, "求和不一致"
```

**2. 与理论值对比:**

```python
rel_diff = abs(beta_total - beta_theory) / beta_theory * 100
if rel_diff > 5.0:
    print(f"警告: 相对偏差 {rel_diff:.2f}% 过大")
```

**3. 网格一致性:**

```python
# 确保通量和共轭通量使用相同网格
with h5py.File('flux_mesh.h5', 'r') as f1, \
     h5py.File('adjoint_flux.h5', 'r') as f2:
    shape1 = f1['grid_shape'][:]
    shape2 = f2['grid_shape'][:]
    assert np.array_equal(shape1, shape2), "网格不匹配"
```

### 物理意义验证

**能量守恒:**

```
Σ (β_i × ν_total) ≈ Σ ν_d,i
```

**空间加权正确性:**

```
β_eff 应该接近 (但不完全等于) β_kinetic = Σν_d / ν_total
```

差异来自共轭通量的空间加权效应。

---

## 故障排除

### 常见问题

#### 1. 文件未生成

**症状:**

```
Error: Unable to open flux_mesh.h5
```

**原因:**

- `clutch_on` 未开启
- 模拟未运行到 `n_batches`
- 通量网格初始化失败

**解决:**

```bash
# 检查 settings.xml
grep -A2 "<clutch>" settings.xml

# 检查控制台输出
grep "flux mesh initialized" output.log

# 手动检查文件
ls -lh flux_mesh.h5 adjoint_flux.h5
```

#### 2. 相对偏差过大

**症状:**

```
Relative difference = 15.3%
```

**可能原因:**

- 网格分辨率不足 (尝试增加网格点数)
- 共轭通量收敛不足 (增加 batch 数)
- 材料成分与 U-235 差异大 (Phase 2 解决)

**诊断:**

```python
# 检查通量场质量
with h5py.File('flux_mesh.h5', 'r') as f:
    flux = f['flux_mean'][:]
    print(f"非零元素: {np.count_nonzero(flux)}")
    print(f"最大值: {flux.max():.3e}")
    print(f"平均值: {flux.mean():.3e}")
```

#### 3. 编译错误

**症状:**

```
error: 'BetaEffective' was not declared in this scope
```

**解决:**

```bash
# 确保头文件包含
grep "beta_effective.h" src/simulation.cpp

# 确保 CMakeLists.txt 包含源文件
grep "beta_effective.cpp" CMakeLists.txt

# 重新配置并编译
cd build
cmake ..
make -j24
```

#### 4. HDF5 结构错误

**症状:**

```
KeyError: 'metadata/flux_file'
```

**原因:**
使用了旧版本生成的 `beta_eff.h5`

**解决:**

```bash
# 删除旧文件重新生成
rm beta_eff.h5
./openmc  # 重新运行
```

### 调试技巧

**1. 详细输出模式:**

修改 `beta_effective.cpp` 临时添加:

```cpp
std::cout << "DEBUG: flux non-zero cells = " << flux_data.size() << "\n";
std::cout << "DEBUG: adjoint non-zero cells = " << adjoint_data.size() << "\n";
std::cout << "DEBUG: overlapping cells = " << overlap_count << "\n";
```

**2. 导出中间结果:**

```python
# 检查分子分母
with h5py.File('beta_eff.h5', 'r') as f:
    num = f['diagnostics/numerator'][:]
    den = f['diagnostics/denominator'][0]
    print(f"分子: {num}")
    print(f"分母: {den}")
    print(f"比值: {num / den}")
```

**3. 可视化通量场:**

```python
import matplotlib.pyplot as plt

with h5py.File('flux_mesh.h5', 'r') as f:
    indices = f['cell_indices'][:]
    flux = f['flux_mean'][:]
    
# 重构稀疏数据到密集网格
dense_flux = np.zeros((19, 19, 19))
for idx, val in zip(indices, flux):
    dense_flux.flat[idx] = val

# 绘制中心切片
plt.imshow(dense_flux[:, :, 9], cmap='hot')
plt.colorbar(label='Flux')
plt.title('Flux at z=0 plane')
plt.show()
```

---

## 参考文献

1. **Duderstadt & Hamilton** (1976), *Nuclear Reactor Analysis*, Wiley.
   - 第6章: 时空动力学，β_eff 的微扰理论推导

2. **Stacey, W.M.** (2007), *Nuclear Reactor Physics*, 2nd Ed.
   - 第9.3节: 有效缓发中子份额的重要性函数定义

3. **OECD/NEA** (2008), *Uncertainty and Target Accuracy Assessment for Innovative Systems Using Recent Covariance Data Evaluations*.
   - β_eff 不确定度传播方法

4. **OpenMC Documentation**: <https://docs.openmc.org>
   - Clutch 方法原理

---

## 贡献者

- **Yaoyang Ren** - 初始实现 (Phase 1)
- 基于 OpenMC 开源项目

## 更新日志

### v1.0 (2025-11-26)

- ✅ Phase 1: 固定 U-235 核数据单能群计算
- ✅ HDF5 结构化输出
- ✅ 理论值验证
- ✅ 稀疏存储优化

### 未来计划

- 🔄 Phase 2: 材料依赖核数据 (预计 2026 Q1)
- 🔄 Phase 3: 多能群扩展 (预计 2026 Q2)
- 🔄 Phase 4: 统计不确定度 (预计 2026 Q3)

---

**License**: MIT (遵循 OpenMC 项目许可证)
