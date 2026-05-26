# Birth-energy beta_eff 方法实现计划

日期：2026-05-26

## 目标

新增一条以裂变中子出生能量为能群语义的 `beta_eff` 方法，并将其作为主输出。
现有 CLUTCH family/group Method E 保留为 `method_comparison/method_e`，用于对比和诊断。

新方法使用 `fission_matrix.h5/adjoint_source_grouped` 中的
`I*(cell,g_birth)`，在 `BetaEffective` 阶段显式折叠
`chi_prompt(g_birth)` 和 `chi_delayed_k(g_birth)`。

## 主公式

```text
D = sum_c V * sum_gin sum_gb phi(c,gin) * Sigma_f(c,gin)
    * [nu_prompt(c,gin) * chi_prompt(c,gb) * I*(c,gb)
       + sum_k nu_delayed_k(c,gin) * chi_delayed_k(c,gb) * I*(c,gb)]

N_k = sum_c V * sum_gin sum_gb phi(c,gin) * Sigma_f(c,gin)
      * nu_delayed_k(c,gin) * chi_delayed_k(c,gb) * I*(c,gb)

beta_k = N_k / D
```

其中：

- `gin` 是诱发裂变的入射能群。
- `gb` 是裂变中子出生能群。
- `I*(cell,gb)` 来自 `fission_matrix.h5/adjoint_source_grouped`。

## 关键改动

1. `BetaEffective::compute_from_files` 增加可选输入
   `fission_matrix_file = "fission_matrix.h5"`。
2. 新增读取逻辑，读取并校验：
   - `shape`
   - `pitch` attribute
   - `n_source_groups` attribute
   - `source_energy_edges`
   - `adjoint_source_grouped`
3. 新增 birth 方法缓存：
   - `birth_adjoint_source_grouped_[cell * G + g_birth]`
   - `birth_source_energy_edges_`
   - `birth_source_n_groups_`
4. 新增计算函数：
   - `compute_denominator_birth_spectrum(...)`
   - `compute_delayed_numerator_birth_spectrum(k, ...)`
5. 顶层 `beta_i`、`beta_total`、`diagnostics/numerator` 和
   `diagnostics/denominator` 切换为 birth-energy 方法结果。
6. 当前 Method E 保留在 `method_comparison/method_e/*`。

## 输出和元数据

主输出保持：

```text
/beta_i
/beta_total
/diagnostics/numerator
/diagnostics/denominator
```

新增或修改 metadata：

```text
metadata/method = birth_spectrum_source_importance
metadata/importance_quantity = fission-matrix source-state importance I*(cell,g_birth)
metadata/group_meaning = fission neutron birth energy group
metadata/chi_birth_spectrum_used = true
metadata/delayed_fraction_location = beta_effective_birth_source_operator
metadata/fission_matrix_file = fission_matrix.h5
metadata/batchwise_ratio_uncertainty = false
```

新增对比输出：

```text
method_comparison/birth_spectrum_source_importance/*
method_comparison/method_e/*
```

在 birth 方法下：

```text
diagnostics/numerator_by_group_energy[k,g]
diagnostics/denominator_by_group_energy[g]
```

表示出生能群贡献。

## 校验规则

- `fission_matrix.h5/shape` 必须与 `flux_mesh.h5/grid_shape` 一致。
- `fission_matrix.h5` 的 `pitch` attribute 必须与 `flux_mesh.h5/grid_pitch[0]` 一致。
- `n_source_groups` 必须等于 `flux_mesh.h5/n_groups`。
- `source_energy_edges` 必须与 `flux_mesh.h5/energy_edges` 一致。
- `adjoint_source_grouped` 尺寸必须等于 `nx * ny * nz * n_groups`。
- 对裂变 cell，必须具备 `sigma_f_groups`、`nu_prompt_groups`、
  `nu_delayed_groups`、`chi_prompt_groups`、`chi_delayed_groups` 和
  `flux_group_map_[cell]`。

## 测试计划

1. 编译验证：

```text
cmake --build build --target openmc -j 4
```

2. HDF5 读取测试：
   - 缺失 `adjoint_source_grouped` 报明确错误。
   - `n_source_groups` 与 `flux n_groups` 不一致时报错。
   - `source_energy_edges` 与 `flux energy_edges` 不一致时报错。
   - `adjoint_source_grouped` 尺寸不匹配时报错。
3. 数值一致性测试：
   - `sum(beta_i) == beta_total`。
   - `denominator > 0`。
   - 所有 `numerator[k] >= 0`。
   - 能群诊断求和等于总 numerator/denominator。
4. 物理回归测试：
   - 均匀单材料模型中 `beta_total` 与材料 delayed fraction 同量级。
   - 各 delayed group 非负。
5. 对比测试：
   - 同一算例同时输出 birth 方法和 Method E。
   - metadata 能明确区分 birth-energy group 和 response-collision group。

## 假设

- 新增方法完成后，顶层 `beta_i` / `beta_total` 使用 birth-energy 方法。
- `adjoint_source_grouped[cell,g]` 视为裂变中子出生状态重要性。
- 第一版要求 `gin` 和 `g_birth` 使用同一套多群边界。
- 第一版不实现 batch-wise ratio uncertainty。
- Method E 不删除，继续作为 response-collision-group CLUTCH 近似对比输出。
