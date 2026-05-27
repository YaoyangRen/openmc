# beta_eff 代码计算流程整理

日期：2026-05-25

更新：2026-05-27

当前 `BetaEffective` 主路径是 birth-energy source-state importance 方法。
它读取 `fission_matrix.h5/adjoint_source_grouped` 得到 `I*(cell,g_birth)`，
并在 `BetaEffective` 阶段显式折叠材料 `chi_prompt(g_birth)` 和
`chi_delayed_k(g_birth)`。现有 CLUTCH family/group Method E 保留为
`method_comparison/method_e` 对比输出。

主公式为：

```text
D = sum_cell volume * sum_gin sum_gbirth phi_gin * sigma_f_gin
    * [nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
       + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]

W_t(cell,gin) = sum_gbirth [
      nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
    + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]

N_k = sum_cell volume * sum_gin phi_gin * sigma_f_gin
      * (nu_delayed_k_gin / nu_total_gin) * W_t(cell,gin)

beta_k = N_k / D
```

其中 `gin` 是诱发裂变入射能群，`gbirth` 是裂变中子出生能群。

## 1. 总体链路

`clutch_on` 打开后，主链路为：

```text
Particle transport
  |
  |-- FluxMesh accumulates track-length forward flux
  |-- FissionMatrix estimates source-state importance I*(source_state)
  |-- GreenFunctionMesh records family-resolved transfer functions
  v
inactive end
  |
  |-- FissionMatrix::compute_adjoint_source()
  |-- write fission_matrix.h5
  v
last batch end
  |
  |-- GreenFunctionMesh::finalize_greenfunction_mesh()
  |-- write transfer_function_data.h5
  |
  |-- AdjointFlux::compute_from_files()
  |-- read transfer_function_data.h5 + fission_matrix.h5
  |-- write adjoint_flux.h5
  |
  |-- FluxMesh::finalize()
  |-- write flux_mesh.h5
  |
  |-- BetaEffective::compute_from_files()
  |-- read flux_mesh.h5 + adjoint_flux.h5 + fission_matrix.h5
  |-- write beta_eff.h5
```

## 2. 输入数据

`flux_mesh.h5` 提供正向多群通量：

```text
flux_group_map_[cell][gin] = phi_gin(cell)
```

`fission_matrix.h5` 提供出生状态重要性：

```text
adjoint_source_grouped[cell * n_groups + gbirth] = I*(cell,gbirth)
source_energy_edges
n_source_groups
shape
pitch
```

`adjoint_flux.h5` 仍会读取，用于 Method E 对比输出。它的
`family_resolved/*/flux_group_mean` 的 `group` 是响应侧诱发裂变碰撞能群，
不是出生能群。

## 3. BetaEffective 主流程

入口：

```cpp
BetaEffective::compute_from_files(
  "flux_mesh.h5", "adjoint_flux.h5", "beta_eff.h5")
```

步骤：

1. 读取 `flux_mesh.h5` 的标量和多群正向通量。
2. 读取 `adjoint_flux.h5` 的 response-weighted importance，用于 Method E 对比。
3. 校验 `flux_mesh.h5` 和 `adjoint_flux.h5` 的网格、pitch、能群数和能群边界一致。
4. 读取 `fission_matrix.h5/adjoint_source_grouped`，并校验 `shape`、`pitch`、`n_source_groups` 和 `source_energy_edges`。
5. 通过几何采样构建 cell 到材料核数据映射。
6. 先计算 Method E 对比结果并写入 `method_comparison/method_e`。
7. 再计算 birth-energy 主结果并写入顶层 `beta_i` / `beta_total`。
8. 写出 `beta_eff.h5`。

## 4. 主算法

主分母函数：

```cpp
compute_denominator_birth_spectrum(flux, volume)
```

离散公式：

```text
D = sum_cell volume * sum_gin sum_gbirth phi_gin * sigma_f_gin
    * [nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
       + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]

W_t(cell,gin) = sum_gbirth [
      nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
    + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]
```

主分子函数：

```cpp
compute_delayed_numerator_birth_spectrum(k, flux, volume)
```

离散公式：

```text
N_k = sum_cell volume * sum_gin phi_gin * sigma_f_gin
      * (nu_delayed_k_gin / nu_total_gin) * W_t(cell,gin)
```

## 5. Method E 对比路径

Method E 仍按响应侧碰撞能群计算：

```text
D_e   = sum_cell volume * sum_g I_total(cell,g)
                            * nu_total_g * sigma_f_g * phi_g

N_e,k = sum_cell volume * sum_g I_delayed_k(cell,g)
                            * nu_total_g * sigma_f_g * phi_g
```

这里 `I_delayed_k(cell,g)` 已经在 transfer-function family 轴包含
`nu_delayed_k / nu_total` 分份，所以 Method E 后处理仍不再乘 `chi_d,k`。

## 6. 输出

`beta_eff.h5` 主结果：

```text
beta_i
beta_total
uncertainty
diagnostics/numerator
diagnostics/denominator
diagnostics/numerator_by_group_energy[k,g_birth]
diagnostics/denominator_by_group_energy[g_birth]
```

关键 metadata：

```text
metadata/method = birth_spectrum_source_importance
metadata/fission_matrix_file = fission_matrix.h5
metadata/importance_quantity = fission-matrix source-state importance I*(cell,g_birth)
metadata/group_meaning = fission neutron birth energy group for I*(cell,g_birth)
metadata/chi_birth_spectrum_used = true
metadata/delayed_fraction_location = yield_fraction_after_birth_importance_folding
metadata/delayed_chi_separate_importance_weight_used = false
metadata/batchwise_ratio_uncertainty = false
```

对比输出：

```text
method_comparison/birth_spectrum_source_importance/*
method_comparison/method_e/*
```

## 7. 当前重要事实

1. 顶层 `beta_i` / `beta_total` 使用 birth-energy 方法。
2. `method_comparison/method_e` 保留 response-collision-group Method E。
3. `diagnostics/*_by_group_energy` 在主结果下表示出生能群贡献。
4. 第一版要求 `gin` 和 `g_birth` 使用同一套多群边界。
5. `uncertainty` 仍为占位零值，尚未实现 batch-wise ratio uncertainty。
