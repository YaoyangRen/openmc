# beta_eff 代码计算流程整理

日期：2026-05-25

更新：2026-05-26

当前 `BetaEffective` 主路径是 CLUTCH family/group Method E。关键语义是：

- `adjoint_flux.h5/family_resolved/*/flux_group_mean` 的 `group` 是响应侧诱发裂变碰撞能群。
- 该 `group` 不是裂变中子出生能群，所以 `BetaEffective` 不引入 `chi_d,k(E_birth)`。
- delayed group 选择性已在 `GreenFunctionMesh` 的 family 轴通过 `nu_d,k / nu_t` 分份进入。
- 因此分子和分母都乘同一个总裂变产生项 `nu_total_g * Sigma_f_g * phi_g`。

主公式为：

```text
D   = sum_cell volume * sum_g I_total(cell,g)
                            * nu_total_g(cell) * sigma_f_g(cell) * phi_g(cell)

N_k = sum_cell volume * sum_g I_delayed_k(cell,g)
                            * nu_total_g(cell) * sigma_f_g(cell) * phi_g(cell)

I_total(cell,g) = I_prompt(cell,g) + sum_k I_delayed_k(cell,g)
beta_k = N_k / D
```

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
  |-- read flux_mesh.h5 + adjoint_flux.h5
  |-- write beta_eff.h5
```

`BetaEffective` 不直接读取 `fission_matrix.h5`。`fission_matrix.h5` 是生成
`adjoint_flux.h5` 的上游输入。

## 2. FluxMesh 输入

`flux_mesh.h5` 提供：

```text
grid_shape
grid_pitch
grid_lower_left
n_groups
energy_edges
cell_indices
flux_mean
flux_group_mean
```

`BetaEffective::read_flux_data()` 将 `flux_group_mean` 整理为：

```text
flux_group_map_[cell][g] = phi_g(cell)
```

当前主路径强制要求多群正向通量数据。

## 3. AdjointFlux 输入

`AdjointFlux` 读取：

```text
transfer_function_data.h5
fission_matrix.h5
```

并计算 response-weighted importance：

```text
I_response(response_cell, family, group)
  = sum_source_state T(source_state -> response_cell, family, group)
                     * I*(source_state)
```

这里的 `group` 来自响应侧裂变碰撞能量。`family` 为：

```text
0      prompt
1..8   delayed_1..delayed_8
```

`adjoint_flux.h5` 输出：

```text
family_resolved/
  prompt/
    cell_indices
    flux_mean
    flux_group_mean
  delayed_1/
    cell_indices
    flux_mean
    flux_group_mean
  ...
  delayed_8/
    cell_indices
    flux_mean
    flux_group_mean
```

`BetaEffective` 正常路径读取 `family_resolved/*/flux_group_mean`。`flux_mean`
只作为兼容和诊断标量保留；多群输入缺失 `flux_group_mean` 时直接报错。

## 4. GreenFunction family 分份

在 fission tally 中，`GreenFunctionMesh` 对响应侧裂变贡献按 family 分份：

```text
prompt    : contribution * (nu_prompt / nu_total)
delayed_k : contribution * (nu_delayed_k / nu_total)
```

因此 `I_delayed_k(cell,g)` 已包含 delayed group `k` 的产额分份。
`BetaEffective` 后处理阶段不能再乘 `nu_delayed_k` 或 `F_delayed_k`。

## 5. BetaEffective 主流程

入口：

```cpp
BetaEffective::compute_from_files(
  "flux_mesh.h5", "adjoint_flux.h5", "beta_eff.h5")
```

步骤：

1. 读取 `flux_mesh.h5` 的标量和多群正向通量。
2. 读取 `adjoint_flux.h5` 的标量和 family-resolved 多群 importance。
3. 校验网格、pitch、能群数和 `energy_edges` 一致。
4. 通过几何采样构建 cell 到材料核数据映射。
5. 使用 Method E family/group 公式计算 `D`、`N_k`、`beta_k`。
6. 写出 `beta_eff.h5`。

## 6. 主算法

分母函数：

```cpp
compute_denominator_upstream_family(flux, volume)
```

离散公式：

```text
D = sum_cell volume * sum_g I_total(cell,g)
                         * nu_total_g * sigma_f_g * phi_g
```

分子函数：

```cpp
compute_delayed_numerator_upstream_family(k, flux, volume)
```

离散公式：

```text
N_k = sum_cell volume * sum_g I_delayed_k(cell,g)
                           * nu_total_g * sigma_f_g * phi_g
```

这里 `g` 是响应侧诱发裂变碰撞能群，不是出生能群。`chi_prompt` 和
`chi_delayed_k` 不进入这个后处理公式。

## 7. 输出

`beta_eff.h5` 主结果保持：

```text
beta_i
beta_total
uncertainty
diagnostics/numerator
diagnostics/denominator
```

新增或明确的诊断和元数据：

```text
diagnostics/numerator_by_group_energy[k,g]
diagnostics/denominator_by_group_energy[g]
diagnostics/numerator_by_group_energy_sum
diagnostics/denominator_by_group_energy_sum
metadata/method = clutch_family_group_method_e
metadata/beta_eff_formula
metadata/family_group_data_used = true
metadata/delayed_fraction_location = transfer_function_family_axis
metadata/group_meaning = response-side induced-fission collision energy group
metadata/chi_birth_spectrum_used = false
metadata/batchwise_ratio_uncertainty = false
method_comparison/method_e/*
```

## 8. 当前重要事实

1. `BetaEffective` 只直接读取 `flux_mesh.h5` 和 `adjoint_flux.h5`。
2. `family_resolved/*/flux_group_mean` 是主路径必需输入。
3. `family_resolved/*/flux_mean` 不再作为正常多群计算的静默降级路径。
4. `nu_delayed_k` 已在 transfer-function family 轴进入，后处理不重复乘。
5. `chi_d,k(E_birth)` 不进入当前 Method E，因为当前 group 不是 birth energy。
6. `uncertainty` 当前仍为占位零值，尚未实现 batch-wise ratio uncertainty。
