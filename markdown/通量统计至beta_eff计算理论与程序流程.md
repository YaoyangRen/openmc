# 通量统计至 beta_eff 计算的理论基础与程序流程

## 摘要

本文档系统说明当前 OpenMC 分支中从正向通量统计、裂变矩阵构造、传递函数统计、response-weighted importance 计算，直至有效缓发中子份额 `beta_eff` 求解的理论基础、离散化方法、程序实现路径和数据文件接口。该流程的核心思想是将正向输运产生的空间和能量依赖信息投影到统一的动力学网格上，并使用裂变矩阵给出的源状态重要性来构造用于 `beta_eff` 计算的伴随权重。

当前实现的主路径不是离线手动几何处理，而是在 OpenMC 运行过程中复用已经由 XML 输入建立的几何、材料和 universe 数据结构。几何网格由 `settings.xml` 中的 `<kinetics_mesh>` 控制，并统一传递给 `FluxMesh`、`FissionMatrix`、`GreenFunctionMesh`、`AdjointFlux` 和 `BetaEffective`。这保证了各中间 HDF5 文件中的网格索引、几何边界和体积因子具有一致的物理含义。

需要特别指出的是，程序中名为 `adjoint_flux.h5` 的文件实际保存的是由传递函数和裂变矩阵源状态重要性卷积得到的 response-weighted importance field。它不是严格意义上通过伴随输运方程直接求解得到的连续能量伴随通量。主 `beta_eff` 方法当前使用 `fission_matrix.h5` 中的 `adjoint_source_grouped`，即 `I*(cell,g_birth)`，作为裂变中子出生能量状态上的重要性。

## 1. 记号与基本离散化

设动力学网格将计算区域离散为若干体积元：

```text
V_c, c = 0, 1, ..., N_c - 1
```

其中 `c` 是一维 cell index。三维网格索引与一维索引之间的映射为：

```text
c = (ix * ny + iy) * nz + iz
```

该索引约定在 `FluxMesh`、`FissionMatrix`、`GreenFunctionMesh` 和 `BetaEffective` 中保持一致。网格间距为 `h`，体积元为：

```text
Delta V = h^3
```

能量空间离散为 `G` 个能群。本文使用：

```text
g_in      诱发裂变的入射中子能群
g_birth   裂变出生中子能群
k         缓发中子先驱核群，k = 1, ..., 8
s         源状态 source_state = cell * G + g_source
```

材料相关核数据包括：

```text
Sigma_f(c,g)          裂变截面
nu_total(c,g)         总裂变中子产额
nu_prompt(c,g)        瞬发中子产额
nu_delayed,k(c,g)     第 k 组缓发中子产额
chi_prompt(c,g_birth) 瞬发裂变谱
chi_delayed,k(c,g_birth) 第 k 组缓发裂变谱
```

所有 `beta_eff` 主公式中的通量均采用体积归一化后的网格通量 `phi(c,g)`。

## 2. 动力学共享网格

### 2.1 网格来源

动力学网格由 `SharedMeshGrid` 创建。输入由 `settings.xml` 的 `<kinetics_mesh>` 给定：

```xml
<kinetics_mesh>
  <pitch>1.0</pitch>
  <auto_bounds>true</auto_bounds>
  <lower_left>-10 -10 -10</lower_left>
  <upper_right>10 10 10</upper_right>
</kinetics_mesh>
```

参数含义如下：

- `pitch`：网格间距，单位为 cm。
- `auto_bounds`：是否从根 universe 的 bounding box 自动获取几何边界。
- `lower_left` 和 `upper_right`：手动指定边界。当 `auto_bounds=false` 或自动边界无界时使用。

网格边界来源在终端输出中以 `bounds_source` 显示：

```text
XML geometry bounding box (read from geometry/model XML)
manual bounds from settings.xml (auto_bounds=false)
manual bounds from settings.xml (geometry XML bounds infinite)
```

因此，输出中可以明确判断当前网格范围来自 XML 几何还是手动输入。

### 2.2 网格一致性

以下 HDF5 文件均写出同一套网格元数据：

```text
grid_shape
grid_lower_left
grid_upper_right
grid_pitch
```

相关文件包括：

```text
flux_mesh.h5
fission_matrix.h5
transfer_function_data.h5
adjoint_flux.h5
beta_eff.h5
```

`BetaEffective` 读取数据时以 `flux_mesh.h5` 为参考网格，逐项校验其他文件的 shape、lower left、upper right 和 pitch。若网格不一致，程序立即报错。这一设计避免了不同物理量在空间索引上错位。

## 3. 正向通量统计

### 3.1 理论定义

连续输运理论中，标量通量可理解为单位体积内径迹长度的期望。对网格体积元 `V_c` 和能群 `g`，track-length estimator 的离散形式可写为：

```text
phi_b(c,g) = 1 / Delta V * sum_{tracks in batch b} w_p * l_{p,c} * 1(g_p = g)
```

其中：

- `w_p` 是粒子权重。
- `l_{p,c}` 是粒子在 cell `c` 内的径迹长度。
- `1(g_p = g)` 表示粒子在该径迹段上的能群归属。
- `b` 表示 batch。

多 batch 统计量为：

```text
phi(c,g) = 1 / N_b * sum_b phi_b(c,g)
```

标量通量为：

```text
phi(c) = sum_g phi(c,g)
```

### 3.2 程序实现

正向通量由 `FluxMesh` 统计。粒子在 `Particle::event_advance()` 中移动前保存起点、方向、能量和能群，然后调用：

```text
FluxMesh::accumulate_track(start_position, direction, weight, distance, energy, group)
```

`accumulate_track()` 的核心步骤为：

1. 将飞行段 `[0, distance]` 裁剪到网格包围盒内。
2. 从裁剪后的起点开始，逐个计算与 x、y、z 三个方向下一网格面的交点。
3. 将跨越多个网格单元的飞行段切分为若干子段。
4. 对每个子段，将 `weight * segment_length` 累积到对应 cell 和能群。

这样避免了将一整段跨网格飞行错误地记入终点 cell 的问题。

### 3.3 批统计和 HDF5 输出

每个 batch 结束时：

```text
FluxMesh::end_batch()
```

会合并线程局部通量缓存，并将本 batch 的 cell-wise track length 累加到跨 batch 统计量中。最终：

```text
FluxMesh::finalize()
```

写出 `flux_mesh.h5`，主要字段为：

```text
cell_indices
flux_mean
flux_std
flux_group_mean
flux_group_std
flux_mean_dense
flux_group_mean_dense
n_groups
energy_edges
grid_shape
grid_lower_left
grid_upper_right
grid_pitch
```

其中 `flux_mean` 和 `flux_group_mean` 已除以 cell volume，因此可直接作为 `beta_eff` 公式中的 `phi(c)` 和 `phi(c,g)` 使用。

当前 `flux_std` 的实现是基于 batch 值的标准差形式，不是严格的均值标准误。若后续需要统计不确定度，应进一步输出 `std / sqrt(N_b)` 或 batch-wise ratio estimator。

## 4. 裂变矩阵与源状态重要性

### 4.1 源状态定义

裂变矩阵以源状态为行索引。源状态定义为：

```text
s = c_source * G + g_source
```

其中 `c_source` 是源粒子出生所在的网格 cell，`g_source` 是源粒子出生能群。源粒子出生时，程序调用：

```text
FissionMatrix::record_source_birth()
```

记录 `source_particle_id -> source_state` 的映射。

### 4.2 裂变事件统计

在裂变碰撞处，程序计算期望裂变中子数：

```text
nu_t = w / k_eff * w_ufs * nu_fission / Sigma_t
```

其中 `w_ufs` 是 uniform fission source weighting 的权重因子。随后调用：

```text
FissionMatrix::record_fission_event()
```

将该贡献累积到：

```text
F(s,j)
```

其中 `s` 为该历史的源状态，`j` 为当前裂变位置所在网格 cell。

对每个源状态进行源粒子数归一化，得到 normalized fission matrix：

```text
M_norm(s,j) = F(s,j) / N_source(s)
```

该矩阵可视为从源状态 `s` 出发，在一个代际传播中对裂变 cell `j` 产生的期望裂变源响应。

### 4.3 经验裂变谱

为了将裂变 cell `j` 上的响应重新投影到出生源状态 `(j,g)`，程序使用源计数构造经验裂变谱：

```text
chi_emp(j,g) = N_source(j,g) / sum_g N_source(j,g)
```

若某个 cell 没有源计数，则采用均匀能群分布作为后备。

### 4.4 伴随源迭代

当前实现将裂变矩阵伴随源解释为源状态重要性 `I*(s)`。伴随迭代采用转置算子形式：

```text
R(j) = sum_g chi_emp(j,g) * I*(j,g)

I_new*(s) = 1 / k_eff * sum_j M_norm(s,j) * R(j)
```

然后对 `I_new*` 做归一化，并以最大范数差值判断收敛：

```text
max_s |I_new*(s) - I_old*(s)| < tolerance
```

该形式对应左特征向量意义上的伴随重要性。它与前向传播形式不同。前向传播会将 `I*(s)` 经 `M_norm(s,j)` 推到响应 cell `j`，而当前伴随形式则用 `M_norm` 的转置作用将响应重要性回传到源状态。

最终输出：

```text
adjoint_source_grouped[s] = I*(cell,g_birth)
adjoint_source[cell] = sum_g I*(cell,g)
```

`beta_eff` 主方法使用的是 `adjoint_source_grouped`。

## 5. 传递函数与 response-weighted importance

### 5.1 传递函数的物理意义

`GreenFunctionMesh` 统计从源状态 `s` 到响应位置 `r` 的裂变响应。对源状态 `s` 和响应 cell `j`，传递函数可表示为：

```text
T_f(s -> j,g)
```

其中：

- `f = 0` 表示 prompt family。
- `f = 1..8` 表示 delayed family 1 到 8。
- `g` 是响应位置处发生裂变碰撞的能群。

程序在 `SCORE_GREENFUNCTION` 中按裂变 family 分配贡献。对于一次裂变碰撞，先计算总期望裂变中子数 `nu_t`，再按产额份额分配：

```text
prompt contribution  = nu_t * nu_prompt / nu_total
delayed k contribution = nu_t * nu_delayed,k / nu_total
```

这意味着 family-resolved 传递函数已经包含了 prompt 或 delayed 产额份额。

### 5.2 源状态归一化

为了使传递函数表示“单位源粒子条件响应”，程序记录每个源状态的源粒子数：

```text
source_state_counts[s]
```

输出 HDF5 时，对每个源状态的累计响应除以该源状态源粒子数：

```text
T_f(s -> j,g) = total_response_f(s -> j,g) / source_state_counts[s]
```

若某源状态存在响应但没有源计数，程序直接报错。这保证了传递函数不是总响应，而是条件平均响应。

### 5.3 transfer_function_data.h5

主要输出包括：

```text
transfer_functions/source_state_*
source_state_indices
source_state_cell
source_state_group
source_counts_per_cell
source_counts_per_state
n_groups
energy_edges
n_families
grid_shape
grid_lower_left
grid_upper_right
grid_pitch
```

每个 `source_state_*` 下保存：

```text
indices
values
source_cell
source_group
source_count
normalized_per_source_particle
```

### 5.4 response-weighted importance

`AdjointFlux` 从 `transfer_function_data.h5` 读取传递函数，从 `fission_matrix.h5` 读取源状态重要性，并计算：

```text
I_response(j,g) = sum_s T(s -> j,g) * I*(s)
```

若存在 family 轴，则同时计算：

```text
I_response,f(j,g) = sum_s T_f(s -> j,g) * I*(s)
```

输出文件为 `adjoint_flux.h5`。需要注意：

```text
I_response(j,g)
```

是 response-weighted importance，不是严格的 transport adjoint flux。程序也在 HDF5 metadata 中写明了这一语义。

`adjoint_flux.h5` 的主要字段为：

```text
cell_indices
flux_mean
flux_group_mean
adjoint_flux_dense
adjoint_flux_group_dense
family_resolved/prompt
family_resolved/delayed_1 ... delayed_8
semantic_metadata
grid_shape
grid_lower_left
grid_upper_right
grid_pitch
```

其中 `family_resolved` 数据主要用于 `beta_eff` 的 Method E 对比路径。

## 6. beta_eff 理论定义

有效缓发中子份额可形式化写为：

```text
beta_eff,k = <F_d,k phi, psi*> / <F phi, psi*>
beta_eff = sum_k beta_eff,k
```

其中：

- `F` 为总裂变源算子。
- `F_d,k` 为第 `k` 组缓发裂变源算子。
- `phi` 为正向通量。
- `psi*` 为伴随重要性函数。

在当前程序中，严格连续形式被离散为网格 cell、入射能群和出生能群上的求和。主方法采用的是裂变矩阵源状态重要性：

```text
I*(c,g_birth)
```

这代表在 cell `c`、出生能群 `g_birth` 上产生一个裂变源中子的相对重要性。

## 7. beta_eff 主方法：birth-energy source-state importance

### 7.1 分母

主方法显式区分诱发裂变能量 `g_in` 和出生能量 `g_birth`。分母为总裂变出生源按重要性加权后的响应：

```text
D =
sum_c Delta V
sum_gin
sum_gbirth
phi(c,g_in) Sigma_f(c,g_in)
[
  nu_prompt(c,g_in) chi_prompt(c,g_birth) I*(c,g_birth)
  +
  sum_k nu_delayed,k(c,g_in) chi_delayed,k(c,g_birth) I*(c,g_birth)
]
```

为简化表达，定义总出生源重要性权重：

```text
W_t(c,g_in) =
sum_gbirth
[
  nu_prompt(c,g_in) chi_prompt(c,g_birth) I*(c,g_birth)
  +
  sum_k nu_delayed,k(c,g_in) chi_delayed,k(c,g_birth) I*(c,g_birth)
]
```

则：

```text
D = sum_c Delta V sum_gin phi(c,g_in) Sigma_f(c,g_in) W_t(c,g_in)
```

### 7.2 分子

第 `k` 个缓发组的分子为：

```text
N_k =
sum_c Delta V
sum_gin
phi(c,g_in) Sigma_f(c,g_in)
[
  nu_delayed,k(c,g_in) / nu_total(c,g_in)
]
W_t(c,g_in)
```

因此：

```text
beta_eff,k = N_k / D
beta_eff = sum_k N_k / D
```

当前实现将 delayed group 的分配放在出生谱重要性折叠之后，即先构造总出生源响应 `W_t`，再乘以 `nu_delayed,k / nu_total`。这样做的理由是 `I*(c,g_birth)` 来自裂变矩阵源状态重要性，不是针对每个 delayed family 的独立输运伴随函数。若再用 delayed chi 单独作为重要性权重，可能对 delayed 组进行重复或过度惩罚。

### 7.3 程序实现

主方法对应以下函数：

```text
BetaEffective::compute_denominator_birth_spectrum()
BetaEffective::compute_delayed_numerator_birth_spectrum()
```

程序步骤为：

1. 读取 `flux_mesh.h5` 中的 `flux_group_mean`，得到 `phi(c,g_in)`。
2. 读取 `fission_matrix.h5` 中的 `adjoint_source_grouped`，得到 `I*(c,g_birth)`。
3. 对每个有通量的 cell，查询几何和材料。
4. 对裂变材料 cell 提取或混合材料核数据。
5. 对 `c`、`g_in`、`g_birth` 和 delayed group 求和。
6. 使用 Kahan summation 降低大量正项累加的浮点误差。
7. 写出 `beta_i` 和 `beta_total`。

## 8. Method E：上游 family-resolved importance 对比方法

除主方法外，程序还保留了 Method E 作为对比。Method E 使用 `adjoint_flux.h5` 中的 `family_resolved` 数据。

### 8.1 分母

定义：

```text
I_total(c,g) = I_prompt(c,g) + sum_k I_delayed,k(c,g)
F_total(c,g) = nu_total(c,g) Sigma_f(c,g) phi(c,g)
```

则：

```text
D_E = sum_c Delta V sum_g I_total(c,g) F_total(c,g)
```

### 8.2 分子

由于 `I_delayed,k` 已经在传递函数中包含了 `nu_delayed,k / nu_total` 的 family 份额，因此分子写为：

```text
N_E,k = sum_c Delta V sum_g I_delayed,k(c,g) F_total(c,g)
```

这里不再额外乘以 delayed yield，以避免重复计数。

### 8.3 适用性

Method E 更接近“上游响应族解析”的解释。它依赖 `transfer_function_data.h5` 中 family 轴的统计质量和 `adjoint_flux.h5` 的 family-resolved group 数据。当前主输出仍采用 birth-energy source-state importance 方法，Method E 作为诊断和对比结果写入 `beta_eff.h5` 的 `method_comparison/method_e`。

## 9. 材料映射与核数据提取

### 9.1 几何查询

`BetaEffective` 不使用手动几何表，而是对每个有通量的网格 cell 执行 OpenMC 几何查询。对 cell 中采样点调用：

```text
exhaustive_find_cell()
```

以获得真实几何 cell 和材料。采样点数量可为：

```text
1, 8, 27
```

单点采样使用 cell center；多点采样使用规则网格点估计 cell 内材料体积分数。

### 9.2 异质 cell 处理

若一个动力学网格 cell 中命中多种裂变材料，程序按采样命中次数计算材料体积分数，并对核数据进行加权混合：

```text
X_cell = sum_m f_m X_m
```

其中 `X` 可以是 `Sigma_f`、`nu_total`、`nu_prompt`、`nu_delayed`、`chi_prompt` 或 `chi_delayed` 等数据。若 cell 只包含非裂变材料、void 或几何查询失败，则不参与 `beta_eff` 的裂变源求和。

### 9.3 核数据

程序通过 `MaterialNuclearDataExtractor` 从 OpenMC 材料和核数据库中提取多群核数据。提取前会根据正向通量谱建立 collapse weights，用于材料数据的能群折合。

## 10. 总体程序流程

完整计算链路可概括如下：

```text
OpenMC XML input
  |
  |-- read settings/materials/geometry/model XML
  |
  |-- SharedMeshGrid::create_from_settings()
  |     |
  |     |-- auto bounds from root universe bounding box
  |     |-- or manual bounds from settings.xml
  |
  |-- transport simulation
        |
        |-- Particle source birth
        |     |-- FissionMatrix::record_source_birth()
        |     |-- GreenFunctionMesh::record_source_birth()
        |
        |-- Particle flight
        |     |-- FluxMesh::accumulate_track()
        |
        |-- Fission collision
        |     |-- FissionMatrix::record_fission_event()
        |     |-- GreenFunctionMesh::accumulate()
        |
        |-- batch end
        |     |-- FluxMesh::end_batch()
        |
        |-- inactive batch end
        |     |-- FissionMatrix::compute_adjoint_source()
        |     |-- write fission_matrix.h5
        |
        |-- final batch end
              |-- write transfer_function_data.h5
              |-- AdjointFlux::compute_from_files()
              |-- write adjoint_flux.h5
              |-- FluxMesh::finalize()
              |-- write flux_mesh.h5
              |-- BetaEffective::compute_from_files()
              |-- write beta_eff.h5
```

## 11. 输入输出文件关系

### 11.1 flux_mesh.h5

提供正向通量：

```text
phi(c)
phi(c,g)
```

并作为 beta 计算的参考网格。

### 11.2 fission_matrix.h5

提供：

```text
M_norm(source_state, fission_cell)
adjoint_source(cell)
adjoint_source_grouped(cell,g_birth)
source_counts
source_energy_edges
```

其中 `adjoint_source_grouped` 是 beta 主方法的核心重要性数据。

### 11.3 transfer_function_data.h5

提供条件响应：

```text
T(source_state -> response_cell, family, response_group)
```

该响应已按 `source_state_counts` 归一化。

### 11.4 adjoint_flux.h5

提供 response-weighted importance：

```text
I_response(c,g)
I_response,family(c,g)
```

其中 family-resolved 数据用于 Method E。

### 11.5 beta_eff.h5

最终输出：

```text
beta_i
beta_total
diagnostics/numerator
diagnostics/denominator
diagnostics/numerator_by_group_energy
diagnostics/denominator_by_group_energy
method_comparison/birth_spectrum_source_importance
method_comparison/method_e
metadata
```

## 12. 数值一致性与物理合理性

当前流程在以下方面具有一致性：

1. 空间网格一致：所有模块使用同一个 `SharedMeshGrid`。
2. 索引一致：一维索引均采用 `(ix * ny + iy) * nz + iz`。
3. 体积归一化一致：通量输出已除以 `Delta V`，beta 公式中再乘以 `Delta V` 完成空间积分。
4. 源状态一致：裂变矩阵和传递函数均使用 `source_state = cell * G + g_source`。
5. 网格元数据一致：beta 读取阶段统一校验所有 HDF5 文件。
6. 几何现实性一致：材料映射通过 OpenMC 几何查询完成，而不是使用手动输入的材料区间。

但仍需注意以下边界：

1. `adjoint_flux.h5` 的名称沿用历史命名，其物理量是 response-weighted importance。
2. 当前通量不确定度不是严格的 beta_eff 不确定度。
3. Method E 的准确性依赖 family-resolved transfer function 的统计质量。
4. 若几何 bounding box 无限，必须提供显式动力学网格边界。
5. 对强异质几何，`1/8/27` 点采样只是材料体积分数近似，网格尺寸和采样点数会影响材料混合精度。

## 13. 建议的验证项目

建议对该链路进行以下验证：

1. 网格一致性验证：检查 `flux_mesh.h5`、`fission_matrix.h5`、`transfer_function_data.h5`、`adjoint_flux.h5` 的四组网格元数据完全一致。
2. 通量统计验证：构造粒子跨多个网格 cell 的简单算例，确认径迹长度被分段计入。
3. 非立方网格索引验证：使用 `nx != ny != nz` 的网格，确认同一空间点在所有模块中映射到相同 cell index。
4. 源状态归一化验证：检查 `source_counts_per_state` 与各 `source_state_*` 的响应归一化关系。
5. 裂变矩阵伴随验证：对小矩阵人工构造 `M_norm`，比较转置 power iteration 与解析左特征向量。
6. beta_eff 闭合验证：检查 `diagnostics/numerator_by_group_energy_sum` 与 `diagnostics/numerator` 的差异仅为浮点求和误差。
7. 材料映射验证：对已知材料分布的几何，检查 `material_cell_counts`、fissionable cell 数和异质 cell 数是否符合预期。

## 14. 结论

当前实现形成了一条从正向通量到 `beta_eff` 的闭合计算路径。正向通量采用网格径迹长度估计；裂变矩阵提供源状态级别的重要性；传递函数提供从源状态到响应位置的条件响应；response-weighted importance 作为诊断和 Method E 输入；最终 `beta_eff` 主结果由 birth-energy source-state importance 方法给出。

从理论结构看，该方法将 `beta_eff` 的经典双线性泛函形式离散到 cell、入射能群和出生能群上，并通过裂变矩阵伴随源近似伴随重要性。其主要优点是与 OpenMC 的 XML 几何、材料和输运流程保持一致，避免手动几何输入导致的空间错配；其主要限制是 response-weighted importance 并非严格伴随输运解，统计不确定度和强异质材料近似仍需后续进一步量化。
