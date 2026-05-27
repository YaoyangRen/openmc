# beta_eff 理论流程与当前代码流程差异分析

日期：2026-05-25

更新：2026-05-27

当前代码已将 `BetaEffective` 主结果切换为 birth-energy source-state
importance 方法。它使用 `fission_matrix.h5/adjoint_source_grouped` 中的
`I*(cell,g_birth)`，并在后处理阶段显式折叠 `chi_prompt` 和
`chi_delayed_k`。

## 1. 理论定义

理论上，第 `k` 个 delayed group 的有效份额可写为：

```math
\beta_{k,\mathrm{eff}}
=
\frac{\langle \phi^\dagger, F_{d,k}\phi \rangle}
       {\langle \phi^\dagger, F_t\phi \rangle}
```

展开裂变源算符后至少包含两个能量变量：

```text
E  = 诱发裂变的入射能量
E' = 裂变中子出生能量
```

严格形式中 delayed 分子包含：

```text
nu_delayed_k(E) * sigma_f(E) * phi(E)
* chi_delayed_k(E') * phi_dagger(E')
```

## 2. 当前主公式

当前代码用离散 birth-energy source-state importance 近似上述结构：

```text
D = sum_cell volume * sum_gin sum_gbirth phi_gin * sigma_f_gin
    * [nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
       + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]

W_t(cell,gin) = sum_gbirth [
      nu_prompt_gin * chi_prompt_gbirth * I*(cell,gbirth)
    + sum_k nu_delayed_k_gin * chi_delayed_k_gbirth * I*(cell,gbirth)]

N_k = sum_cell volume * sum_gin phi_gin * sigma_f_gin
      * (nu_delayed_k_gin / nu_total_gin) * W_t(cell,gin)
```

其中：

```text
I*(cell,gbirth) = fission_matrix.h5/adjoint_source_grouped[cell,gbirth]
```

因此当前主路径明确区分：

- `gin`：诱发裂变的入射能群。
- `gbirth`：裂变中子出生能群。

## 3. 与 Method E 的关系

`adjoint_flux.h5/family_resolved/*/flux_group_mean` 的 `group` 仍是响应侧诱发裂变碰撞能群，
不是出生能群。因此 Method E 仍按以下公式作为对比路径：

```text
D_e   = sum_cell volume * sum_g I_total(cell,g) * F_total_g
N_e,k = sum_cell volume * sum_g I_delayed_k(cell,g) * F_total_g
```

Method E 中 `I_delayed_k` 已经通过 transfer-function family 轴包含
`nu_delayed_k / nu_total` 分份，所以 Method E 后处理不再乘 `chi_d,k`。

## 4. 已实现内容

- 顶层 `beta_i` / `beta_total` 使用 birth-energy 方法。
- `BetaEffective` 读取并校验 `fission_matrix.h5/adjoint_source_grouped`。
- HDF5 metadata 写入 `method = birth_spectrum_source_importance`。
- HDF5 metadata 写入 `chi_birth_spectrum_used = true`。
- HDF5 metadata 写入 `delayed_fraction_location = yield_fraction_after_birth_importance_folding`。
- HDF5 metadata 写入 `delayed_chi_separate_importance_weight_used = false`。
- 能群诊断 `numerator_by_group_energy` 和 `denominator_by_group_energy` 表示出生能群贡献。
- Method E 保留在 `method_comparison/method_e`。

## 5. 剩余差异

| 项目 | 严格理论流程 | 当前代码流程 | 影响 |
|---|---|---|---|
| 伴随量 | 真伴随通量或等价重要性 | fission-matrix source-state importance | 仍是裂变矩阵近似 |
| 能量变量 | 连续 `E` 和 `E'` | 离散 `gin` 和 `gbirth` | 多群离散误差 |
| 出生谱 | 连续 `chi(E')` | 材料多群 `chi_groups` | 由核数据抽样/折合质量决定 |
| 空间处理 | 连续空间积分 | mesh cell 求和和材料采样 | 有空间离散误差 |
| 不确定度 | ratio estimator | `uncertainty` 仍为零占位 | 尚无统计误差传播 |

## 6. 结论

当前主路径已经从 response-collision-group Method E 切换到
birth-energy source-state importance。现在 `chi_p(E_birth)` 和
`chi_d,k(E_birth)` 先用于折叠总出生源重要性 `W_t`；delayed 组分子再按
`nu_delayed_k / nu_total` 分配同一个 `W_t`。这样避免把
`chi_d,k(E_birth)` 作为单独重要性权重重复压低 delayed 分子；Method E
仍作为响应侧碰撞能群对比结果保留。
