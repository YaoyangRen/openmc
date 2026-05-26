# beta_eff 理论流程与当前代码流程差异分析

日期：2026-05-25

更新：2026-05-26

当前代码实现采用 CLUTCH family/group Method E 作为 `beta_eff.h5` 主结果。
本次修正明确了一个关键差异：`adjoint_flux.h5/family_resolved/*/flux_group_mean`
里的 `group` 是响应侧诱发裂变碰撞能群，不是裂变中子出生能群。因此当前
`BetaEffective` 仍不引入 `chi_d,k(E_birth)`。

## 1. 理论定义

理论上，第 `k` 个 delayed group 的有效份额可写为：

```math
\beta_{k,\mathrm{eff}}
=
\frac{\langle \phi^\dagger, F_{d,k}\phi \rangle}
       {\langle \phi^\dagger, F_t\phi \rangle}
```

若展开连续能量裂变源算符，分子包含两类能量变量：

```text
E  = 诱发裂变的入射能量
E' = 裂变中子出生能量
```

严格后处理形式可以包含：

```text
nu_delayed_k(E) * sigma_f(E) * phi(E) * chi_delayed_k(E') * phi_dagger(E')
```

这要求伴随权重作用在出生中子状态 `E'` 上。

## 2. 当前代码计算对象

当前 `adjoint_flux.h5` 中的 family-resolved 数据来自：

```text
I_response(response_cell, family, group)
  = sum_source_state T(source_state -> response_cell, family, group)
                     * I*(source_state)
```

其中 `group` 是响应侧裂变碰撞能群。它描述的是贡献链到达响应 cell 后，
发生诱发裂变碰撞时的能量分组，不是裂变后新生中子的出生能量分组。

因此当前主公式为：

```text
D   = sum_cell volume * sum_g I_total(cell,g)
                            * nu_total_g * sigma_f_g * phi_g

N_k = sum_cell volume * sum_g I_delayed_k(cell,g)
                            * nu_total_g * sigma_f_g * phi_g
```

其中：

```text
I_total(cell,g) = I_prompt(cell,g) + sum_k I_delayed_k(cell,g)
```

## 3. 为什么不引入 chi_d,k(E_birth)

如果在当前公式中直接乘 `chi_d,k(E_birth)`，会把响应侧碰撞能群 `g` 当作出生能群
`E'` 使用，造成能量语义错位。

当前 delayed group 选择性已经在 `GreenFunctionMesh` 中进入：

```text
prompt    : contribution * (nu_prompt / nu_total)
delayed_k : contribution * (nu_delayed_k / nu_total)
```

所以 `I_delayed_k(cell,g)` 已经携带 delayed group `k` 的分份。后处理分子必须乘
`F_total_g = nu_total_g * sigma_f_g * phi_g`，而不是再乘 `nu_delayed_k`、
`F_delayed_k` 或 `chi_delayed_k`。

## 4. 已实现的修正

当前代码已经按上述语义调整：

- `BetaEffective` 主方法为 `clutch_family_group_method_e`。
- 主公式读取 `family_resolved/*/flux_group_mean`。
- 缺失 family group 数据时直接 `fatal_error()`，不静默退回 `flux_mean` 标量。
- 分母使用 `I_prompt(cell,g) + sum_k I_delayed_k(cell,g)`。
- 分子使用对应 `I_delayed_k(cell,g)`。
- 分子和分母都乘 `nu_total_g * sigma_f_g * phi_g`。
- HDF5 元数据记录 `delayed_fraction_location = transfer_function_family_axis`。
- HDF5 元数据记录 `chi_birth_spectrum_used = false`。
- 能群诊断 `numerator_by_group_energy` 和 `denominator_by_group_energy`
  表示响应侧碰撞能群贡献。

## 5. 当前理论差异

| 项目 | 严格理论流程 | 当前代码流程 | 影响 |
|---|---|---|---|
| 伴随量 | 真伴随通量或等价源状态重要性 | response-weighted importance | 不是严格输运伴随求解 |
| 能量变量 | 区分入射能量 `E` 和出生能量 `E'` | 使用响应侧碰撞能群 `g` | 不能直接使用出生谱 `chi(E')` |
| delayed 信息 | 可在 `F_d,k` 和 `chi_d,k` 中进入 | 在 transfer-function family 轴进入 | 后处理不能重复乘 delayed yield |
| 分子 | `phi_dagger * F_d,k * phi` | `I_delayed_k(cell,g) * F_total_g` | 记账方式不同但语义自洽 |
| 分母 | 总裂变源伴随加权 | family importance 求和后乘 `F_total_g` | 依赖 family 分份一致性 |
| 空间处理 | 连续空间积分 | mesh cell 求和和材料采样 | 有空间离散误差 |
| 不确定度 | ratio estimator | `uncertainty` 仍为零占位 | 尚无统计误差传播 |

## 6. 后续计划

短期计划：

1. 保持 Method E 作为当前主方法，不引入 `chi_d,k(E_birth)`。
2. 增加单元级测试，验证 `I_delayed_k` 已含 delayed 分份时分子只乘 `F_total_g`。
3. 增加 HDF5 读取测试，覆盖缺失 `flux_group_mean`、尺寸不匹配和能群不一致。
4. 增加一致性测试，检查 `sum_k beta_i == beta_total`、分母大于零、分子非负。
5. 增加诊断求和测试，确认能群诊断求和等于总 numerator/denominator。

长期如果要引入 `chi_d,k(E_birth)`，不能在当前 response-group `adjoint_flux`
后处理里直接相乘。需要先建立出生状态重要性语义，例如显式的
`I*(cell,g_birth)` 或等价的 birth-energy adjoint/source-importance 算符，再在
源算符层面折叠 `chi_prompt` 和 `chi_delayed_k`。

## 7. 结论

当前实现和理论公式的主要差异不是简单缺少一个 `chi_d,k` 因子，而是能量变量语义不同。
当前 `group` 是响应侧诱发裂变碰撞能量，delayed group 分份已经在 transfer function
family 轴中完成。因此当前 Method E 中保持不引入 `chi_d,k(E_birth)` 是一致的实现选择。

