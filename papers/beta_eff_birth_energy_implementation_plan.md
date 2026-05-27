# Birth-energy beta_eff 方法实现记录

日期：2026-05-27

## 摘要

本计划已落地为当前 `BetaEffective` 主路径。顶层 `beta_i` / `beta_total`
使用 birth-energy source-state importance 方法；原 CLUTCH family/group
Method E 保留为 `method_comparison/method_e`。

## 已实现主公式

```text
D = sum_c V * sum_gin sum_gb phi(c,gin) * Sigma_f(c,gin)
    * [nu_prompt(c,gin) * chi_prompt(c,gb) * I*(c,gb)
       + sum_k nu_delayed_k(c,gin) * chi_delayed_k(c,gb) * I*(c,gb)]

W_t(c,gin) = sum_gb [
      nu_prompt(c,gin) * chi_prompt(c,gb) * I*(c,gb)
    + sum_k nu_delayed_k(c,gin) * chi_delayed_k(c,gb) * I*(c,gb)]

N_k = sum_c V * sum_gin phi(c,gin) * Sigma_f(c,gin)
      * (nu_delayed_k(c,gin) / nu_total(c,gin)) * W_t(c,gin)

beta_k = N_k / D
```

其中 `gin` 是诱发裂变入射能群，`gb` 是裂变中子出生能群。

## 已实现接口和数据流

- `BetaEffective::compute_from_files` 增加可选参数
  `fission_matrix_file = "fission_matrix.h5"`。
- 新增读取 `fission_matrix.h5`：
  - `shape`
  - `pitch` attribute
  - `n_source_groups` attribute
  - `source_energy_edges`
  - `adjoint_source_grouped`
- 新增主计算函数：
  - `compute_denominator_birth_spectrum`
  - `compute_delayed_numerator_birth_spectrum`
- 顶层结果使用 birth 方法。
- `method_comparison/birth_spectrum_source_importance` 写入主方法副本。
- `method_comparison/method_e` 写入响应侧碰撞能群 Method E 对比结果。

## 输出语义

```text
metadata/method = birth_spectrum_source_importance
metadata/importance_quantity = fission-matrix source-state importance I*(cell,g_birth)
metadata/group_meaning = fission neutron birth energy group for I*(cell,g_birth)
metadata/chi_birth_spectrum_used = true
metadata/delayed_fraction_location = yield_fraction_after_birth_importance_folding
metadata/delayed_chi_separate_importance_weight_used = false
metadata/fission_matrix_file = fission_matrix.h5
```

`diagnostics/numerator_by_group_energy` 和
`diagnostics/denominator_by_group_energy` 在主方法下表示出生能群贡献。

## 校验规则

- `fission_matrix.h5/shape` 必须与 `flux_mesh.h5/grid_shape` 一致。
- `fission_matrix.h5` 的 `pitch` attribute 必须与 `flux_mesh.h5/grid_pitch[0]` 一致。
- `n_source_groups` 必须等于 `flux_mesh.h5/n_groups`。
- `source_energy_edges` 必须与 `flux_mesh.h5/energy_edges` 一致。
- `adjoint_source_grouped` 尺寸必须等于 `nx * ny * nz * n_groups`。

## 后续测试缺口

- 还需要补单元级假数据测试，验证 birth 公式中 `W_t * nu_delayed_k / nu_total` 的使用。
- 还需要补 HDF5 读取错误路径测试。
- 还需要补 batch-wise ratio uncertainty；当前 `uncertainty` 仍为零占位。
