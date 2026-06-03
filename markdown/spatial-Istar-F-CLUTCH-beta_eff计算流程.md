# Spatial-I* F-CLUTCH beta_eff 计算流程与算法说明

本文档记录当前代码中 `beta_eff` 的主计算方法。当前实现已经收敛为文献 CLUTCH/F-CLUTCH 的空间裂变源重要性方法：

```text
source state = cell
importance   = I*(cell)
primary      = fclutch_spatial
```

该方法不再使用出生能量维度的伴随源 `I*(cell, g_birth)`，也不再默认生成 `transfer_function_data.h5`、`adjoint_flux.h5` 或 `flux_mesh.h5`。

## 总体流程

当前 `clutch_on` 下的默认 `beta_eff` 流程为：

```text
inactive batches:
  1. 记录每个源粒子的出生 cell
  2. 记录该源粒子产生的裂变源 site 所在 child cell
  3. 累积 cell -> cell 裂变矩阵
  4. inactive 结束后求解空间伴随裂变源 I*(cell)

active batches:
  5. 对每个实际生成的 fission source site 计分
  6. 按 delayed group 累积分子，按所有 fission source site 累积分母
  7. active 结束后按 batch 保存分子和分母

final:
  8. 用 ratio-of-means 计算 beta_i、beta_total 和不确定度
  9. 写出 beta_eff.h5
```

默认输出只依赖：

```text
fission_matrix.h5
beta_eff.h5
```

`flux_mesh.h5` 仅在用户显式启用 `flux_mesh_on` 时作为诊断输出生成。

## 空间裂变矩阵

空间裂变矩阵由 inactive batches 生成。每个源粒子在开始历史时记录其出生位置：

```text
parent_cell = cell(source birth position)
```

当该历史产生 fission source site 时，记录：

```text
child_cell = cell(fission source site position)
```

矩阵原始计数为：

```text
F_raw(parent_cell, child_cell) += site_weight
```

每个 parent cell 的源粒子权重计数为：

```text
C(parent_cell) += source_weight
```

归一化后的空间裂变矩阵为：

```text
F(parent_cell, child_cell) = F_raw(parent_cell, child_cell) / C(parent_cell)
```

HDF5 中的矩阵语义为：

```text
row = parent_cell
col = child_cell
```

对应源码位置：

```text
FissionMatrix::record_source_birth()
FissionMatrix::record_fission_site()
FissionMatrix::start_new_batch()
```

## 伴随源迭代

inactive batches 结束后，代码对空间裂变矩阵的转置算子做幂迭代，求解：

```text
I*(i) = (1 / k_eff) * sum_j F(i, j) * I*(j)
```

其中：

```text
i = parent cell
j = child cell
F(i, j) = normalized spatial fission matrix
```

每次迭代后对 `I*` 做归一化：

```text
sum_i I*(i) = 1
```

收敛判据为：

```text
max_i |I_new*(i) - I_old*(i)| < adjoint_tolerance
```

相关设置来自：

```text
settings::adjoint_initial_guess
settings::adjoint_max_iterations
settings::adjoint_tolerance
```

当前伴随源输出在：

```text
fission_matrix.h5/adjoint_source
```

其含义为：

```text
Spatial CLUTCH adjoint fission source I*(cell)
```

## Active Batch 事件计分

active batches 中不再做后处理矩阵折叠，而是在实际产生 fission source site 时直接计分。

对每个 fission source site：

```text
c = cell(site.r)
w = site.wgt
g_d = site.delayed_group
I = I*(c)
```

分母计分：

```text
D_b += w * I
```

若该 site 是第 `k` 组 delayed neutron：

```text
N_{k,b} += w * I
```

prompt neutron 只进入分母，不进入 delayed group 分子：

```text
delayed_group == 0 -> prompt
```

delayed group 映射为：

```text
delayed_group = 1..8 -> beta_i[0..7]
```

超过 `[0, 8]` 的 delayed group 不计入分子，并记录 warning 计数。

对应源码位置：

```text
BetaEffectiveAccumulator::score_fission_site()
```

CE 和 MG 的 fission site 均调用同一个 spatial-I* 计分接口：

```text
score_fission_site(site.r, site.wgt, site.delayed_group)
```

## Ratio-of-Means 统计

每个 active batch 结束时保存该 batch 的：

```text
D_b
N_{k,b}
```

最终 delayed group beta 为 ratio-of-means：

```text
beta_k = mean_b(N_{k,b}) / mean_b(D_b)
```

总 beta 为：

```text
beta_total = sum_k beta_k
```

这种形式避免了先对每个 batch 做比值再平均带来的偏差。

## 不确定度估计

对每个 delayed group，使用 delta method，并显式考虑分母涨落：

```text
z_{k,b} = N_{k,b} - beta_k * D_b

sigma(beta_k) =
  sqrt( Var(z_k) / (n_batch * mean(D)^2) )
```

总 beta 的不确定度同理：

```text
N_total,b = sum_k N_{k,b}
z_total,b = N_total,b - beta_total * D_b
```

这比只统计分子方差更合理，因为 `beta_eff` 是分子和分母的比值。

## HDF5 输出

### fission_matrix.h5

主要数据集和属性：

```text
row_indices
col_indices
data_raw
data_normalized
source_counts
adjoint_source
```

关键属性：

```text
filetype = fission_matrix_sparse
version = 5.0
matrix_semantics = row=parent_cell, col=child_cell
n_source_groups = 1
adjoint_source_description = Spatial CLUTCH adjoint fission source I*(cell)
```

### beta_eff.h5

根数据集保持兼容：

```text
beta_i
beta_total
uncertainty
```

元数据：

```text
metadata/primary_method = fclutch_spatial
metadata/theory_reference = Qiu2016_F_CLUTCH_Eq31_Eq43
metadata/uses_fission_event_sites = 1
metadata/uses_birth_energy_importance = 0
metadata/source_state_definition = cell
```

方法组：

```text
method/fclutch_spatial/beta_i
method/fclutch_spatial/beta_total
method/fclutch_spatial/numerator
method/fclutch_spatial/denominator
method/fclutch_spatial/uncertainty
method/fclutch_spatial/beta_total_uncertainty
method/fclutch_spatial/batch_ids
method/fclutch_spatial/batch_denominator
method/fclutch_spatial/batch_numerator
```

诊断组：

```text
diagnostics/total_fission_sites
diagnostics/total_scored_sites
diagnostics/total_dropped_sites
diagnostics/total_invalid_delayed_group_sites
diagnostics/active_batches_scored
diagnostics/numerator
diagnostics/denominator
```

## 已移除的默认路径

当前 `clutch_on` 默认路径不再执行：

```text
GreenFunctionMesh creation
transfer_function_data.h5 write
AdjointFlux::compute_from_files()
adjoint_flux.h5 write
birth-energy I*(cell,g_birth) beta_eff
```

也不再因为 CLUTCH tally 自动开启：

```text
settings::flux_mesh_on = true
```

因此默认 `beta_eff` 运行开销主要来自：

```text
inactive fission matrix accumulation
spatial adjoint source iteration
active fission-site beta_eff scoring
```

## 与理论公式的对应

连续形式的有效缓发中子份额可写为：

```text
beta_eff = <psi*, F_delayed phi> / <psi*, F_total phi>
```

当前 F-CLUTCH 事件估计中，`psi*` 由空间裂变源重要性 `I*(cell)` 近似。对实际 fission source site 计分后：

```text
D   = sum_all_fission_sites w_site * I*(cell_site)
N_k = sum_delayed_group_k_sites w_site * I*(cell_site)
```

最终：

```text
beta_k = N_k / D
beta_eff = sum_k beta_k
```

这对应文献中用 fission-event/progenitor-set 权重估计反应率比值的思想；当前实现选择的 source state 是空间 cell，而不是 `(cell, birth energy)`。

## 当前方法的适用性

该方法适合当前问题的原因：

```text
1. cell 状态空间远小于 cell * energy group
2. inactive 裂变矩阵更容易得到非零、连通的伴随源分布
3. 避免 birth-energy 高维伴随源欠采样导致 beta_eff 系统性偏小
4. 已在当前算例中给出接近 MCNP 参考的总 beta_eff
```

需要注意：

```text
1. 分组 beta_i 仍依赖 delayed group 采样和核数据
2. 若 adjoint iteration 未收敛，应增加 adjoint_max_iterations 或适当放宽 tolerance
3. 若 inactive 统计不足，I*(cell) 仍可能稀疏，需要增加 inactive histories 或粗化 kinetics mesh
```
