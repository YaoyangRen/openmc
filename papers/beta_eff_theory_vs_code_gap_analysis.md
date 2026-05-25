# beta_eff 理论流程与当前代码流程差异分析

日期：2026-05-25

本文对比对象：

- 理论流程：从 `beta_eff` 的伴随加权定义和 CLUTCH-like 贡献链接流程出发。
- 代码流程：当前 `flux_mesh.h5 -> adjoint_flux.h5 -> beta_eff.h5` 的 Method E 实现。

这里不评价代码“对错”，只明确两者计算对象、数据语义和近似位置的差异。

## 1. 总体差异

理论上，`beta_eff` 计算的是：

```math
\beta_{k,\mathrm{eff}}
=
\frac{\langle \phi^\dagger, F_{d,k}\phi \rangle}
       {\langle \phi^\dagger, F_t\phi \rangle}
```

也就是缓发第 `k` 群裂变中子源与总裂变中子源的伴随加权贡献比。

当前代码实际主公式是 Method E：

```text
D   = sum_cell volume * I_total(cell)     * F_total(cell)
N_k = sum_cell volume * I_delayed_k(cell) * F_total(cell)
beta_k = N_k / D
```

其中：

```text
I_total(cell) = I_prompt(cell) + sum_k I_delayed_k(cell)
F_total(cell) = sum_g nu_total_g * sigma_f_g * phi_g
```

核心差异是：理论式中 delayed 选择性在 `F_{d,k}` 或 `chi_{d,k}` 中显式出现；当前代码把 delayed 选择性提前放进了上游 `I_delayed_k`，因此后处理阶段统一乘 `F_total`。

## 2. 伴随量语义不同

### 理论流程

理论流程假设有真正的伴随通量或等价的伴随重要性：

```math
\phi^\dagger(r,E)
```

它应该能作为任意裂变中子源项的权重，用来衡量该源中子对未来裂变链或响应的贡献。

如果写成连续能量形式，伴随权重一般作用在裂变中子出生状态上：

```math
\chi_{d,k}(E')\phi^\dagger(r,E')
```

这里 `E'` 是出生中子能量。

### 当前代码流程

当前 `adjoint_flux.h5` 中的量来自：

```text
Phi_dag(response_cell, family, group)
  = sum_source_state T(source_state -> response_cell, family, group)
                     * I*(source_state)
```

代码元数据也明确写了它是：

```text
Response-weighted importance field, NOT true adjoint flux
```

而且 group 的语义是响应位置发生裂变的碰撞能量，不是裂变中子出生能量。

### 差异影响

理论流程里的 `phi^\dagger(r,E')` 和当前代码里的 `Phi_dag(cell, group)` 不是同一个严格物理量。

因此：

- 理论式可以自然使用 `chi_p(E')`、`chi_d,k(E')`；
- 当前代码若直接用 response collision group 去乘出生谱 `chi`，会发生能量语义错位；
- 当前 Method E 避开了这个问题，把 delayed group 选择性放在 `nu_d,k(E_collision) / nu_total(E_collision)` 上。

## 3. delayed group 信息进入位置不同

### 理论流程

理论上，分子可以写成：

```math
N_k =
\int \phi(r,E)\Sigma_f(r,E)\nu_{d,k}(r,E)
\chi_{d,k}(E')\phi^\dagger(r,E')\,dr\,dE\,dE'
```

delayed group 信息主要通过两部分进入：

- `nu_{d,k}(E)`：第 `k` 群 delayed yield；
- `chi_{d,k}(E')`：第 `k` 群 delayed neutron birth spectrum。

### 当前代码流程

当前代码在 `SCORE_GREENFUNCTION` 阶段就按 family 分份：

```text
prompt:    nu_t * (nu_p  / nu_total_yield)
delayed_k: nu_t * (nu_dk / nu_total_yield)
```

也就是说，`I_delayed_k(cell)` 已经包含 delayed-k 的产额分份。

所以 `beta_effective.cpp` 中：

```text
N_k = sum_cell volume * I_delayed_k(cell) * F_total(cell)
```

而不是：

```text
N_k = sum_cell volume * I_delayed_k(cell) * F_delayed_k(cell)
```

### 差异影响

这是两套等价意图下的不同记账方式，但不能混用：

- 理论后处理式：`I*` 不含 delayed group 信息，则分子乘 `F_delayed_k`。
- 当前代码式：`I_delayed_k` 已含 delayed group 分份，则分子乘 `F_total`。

若在当前代码里再乘 `nu_d,k` 或 `F_delayed_k`，会重复计数 delayed yield。

## 4. 能量维度处理不同

### 理论流程

理论流程至少有两个能量变量：

```text
E  = 诱发裂变的入射能量
E' = 裂变中子出生能量
```

严格公式中：

- `nu(E)`、`sigma_f(E)` 依赖入射能量；
- `chi(E')` 和伴随权重依赖出生能量；
- 两者通过裂变源算符连接。

### 当前代码流程

当前代码有两类能量维度：

1. `source_state = source_cell * n_source_groups + g_source`

用于 `FissionMatrix` 和 `GreenFunctionMesh` 的源状态。

2. `group` in transfer function

用于响应位置裂变贡献的能群，来自 `GreenFunctionMesh::accumulate()` 的 `energy_eV` 或 `mg_group`。

`AdjointFlux` 卷积后写出的 group data 语义是：

```text
response collision energy group
```

`BetaEffective` 计算 `F_total` 时使用正向通量分群：

```text
F_total(cell) = sum_g nu_total_g * sigma_f_g * phi_g
```

但当前 Method E 对 `I_delayed_k` 使用的是 `flux_mean` 标量 family importance，没有用 `family_resolved/*/flux_group_mean`。

### 差异影响

理论流程可以区分：

```text
入射能量 E
出生能量 E'
```

当前代码把后处理阶段的 family importance 折叠成 cell 标量：

```text
I_delayed_k(cell)
```

再乘多群 `F_total(cell)`。

这会丢掉 family importance 在响应能群上的结构。更接近当前代码数据语义的多群形式应是：

```text
D = sum_cell volume * sum_g I_total(cell,g)
                          * nu_total_g * sigma_f_g * phi_g
```

```text
N_k = sum_cell volume * sum_g I_delayed_k(cell,g)
                            * nu_total_g * sigma_f_g * phi_g
```

而当前实际是：

```text
D = sum_cell volume * I_total(cell)
                  * sum_g nu_total_g * sigma_f_g * phi_g
```

```text
N_k = sum_cell volume * I_delayed_k(cell)
                    * sum_g nu_total_g * sigma_f_g * phi_g
```

两者的差别是是否保留 `I(cell,g)` 与 `F_total(cell,g)` 的能群相关性。

## 5. 空间积分离散方式不同

### 理论流程

理论式是连续空间积分：

```math
\int_V (...) \, dr
```

材料、通量、伴随量和核数据都应在空间点上定义。

### 当前代码流程

当前代码在规则网格 cell 上计算：

```text
sum_cell volume * ...
```

每个 mesh cell 的材料核数据由 `build_cell_material_map()` 几何采样得到。

默认：

```text
n_sample_points = 27
weighting_mode = REACTION_RATE_WEIGHTED
```

异质 cell 会被折算成一个混合的 `MaterialNuclearData`。

### 差异影响

理论流程是点态积分；当前代码是网格均匀化积分。

误差主要来自：

- mesh cell 内通量变化被平均；
- mesh cell 内材料异质性被采样近似；
- 混合材料核数据不一定等价于严格的空间积分。

这类误差在强异质几何、边界区域、燃料-慢化剂交界处会更明显。

## 6. 统计处理不同

### 理论流程

理论上 `beta_eff` 是 ratio estimator：

```text
beta_k = N_k / D
```

统计不确定度应考虑 `N_k` 与 `D` 的相关性，通常需要 batch-wise 的 numerator/denominator 历史。

### 当前代码流程

当前 `beta_eff.h5` 中：

```text
uncertainty = 8 个 0
```

代码只写最终累积值，没有在 `BetaEffective` 中做 batch-wise ratio uncertainty propagation。

### 差异影响

当前输出没有真实 `beta_i` 或 `beta_total` 的统计不确定度。

如果后续要和 MCNP/Serpent/实验值严肃对比，需要补：

- 每个 batch 的 `N_k`；
- 每个 batch 的 `D`；
- ratio estimator 方差；
- delayed group 之间和 numerator/denominator 的协方差处理。

## 7. 分母定义上的差异

### 理论流程

理论分母是全部裂变中子源的伴随加权贡献：

```math
D = \langle \phi^\dagger, F_t\phi \rangle
```

若严格分 prompt/delayed birth spectrum，分母的伴随权重应与总裂变源算符一致。

### 当前代码流程

当前分母使用：

```text
I_total(cell) = I_prompt(cell) + sum_k I_delayed_k(cell)
D = sum_cell volume * I_total(cell) * F_total(cell)
```

其中 `I_prompt` 和 `I_delayed_k` 来自上游 family-resolved transfer function。

### 差异影响

当前分母不是用一个独立求得的 `phi_dagger_total(cell,g)` 直接进入严格裂变源算符，而是把 family-resolved importance 求和后作为总 importance。

如果 family 划分严格守恒，`I_prompt + sum I_delayed_k` 应接近 total importance；但由于当前后处理读取的是各 family 的 `flux_mean`，不是按 group 与 `F_total_g` 同步相乘，仍可能存在能群折叠误差。

## 8. 源端重要性处理不同

### 理论流程

CLUTCH-like 理论流程需要源状态重要性：

```text
I*(source state)
```

它应表示源中子对未来响应的贡献能力。

### 当前代码流程

当前 `FissionMatrix::compute_adjoint_source()` 先得到空间上的 `q(j)`，再用经验源能谱分配到 group：

```text
I*_new(j,g) = q(j) * chi_empirical(g|j)
```

其中 `chi_empirical(g|j)` 来自源计数：

```text
source_counts[j,g] / sum_g source_counts[j,g]
```

### 差异影响

当前 source-state 重要性不是直接由完整 `(cell,g) -> (cell,g)` 裂变源转移严格闭合出来的，而是在幂迭代中用 `q(j)` 乘经验能谱展开。

这比纯空间 `I*(cell)` 更细，但仍是近似：

- 源能群结构主要由经验源分布给出；
- 能群之间对未来响应的差异未必完全由 `chi_empirical` 捕获；
- 若源能谱统计不足，会影响 `I*(cell,g)`。

## 9. 核数据使用方式不同

### 理论流程

理论流程中，`nu_t`、`nu_d,k`、`sigma_f`、`chi` 都应在对应材料、入射能量和出生能量上连续处理。

### 当前代码流程

当前 `BetaEffective` 从 OpenMC 材料/核素数据中抽取多群核数据，并放入 `MaterialNuclearData`。

计算 Method E 时真正进入主公式的是：

```text
nu_total_groups[g]
sigma_f_groups[g]
flux_group_map[cell][g]
```

`nu_delayed_groups` 在 Method E 后处理公式中不直接使用，因为 delayed 分份已经在上游 transfer function 中进入。

### 差异影响

当前代码把 delayed yield 的使用分成两处：

- 上游 transfer function：用事件核素和事件能量的 `nu_d,k / nu_total`；
- 后处理 `BetaEffective`：只使用 `nu_total_g * sigma_f_g * phi_g`。

这和“后处理直接用 `nu_d,k` 构造分子”的理论写法不同，但与当前 family-resolved 记账方式一致。

## 10. 输出对象差异

### 理论流程

理论上输出应包括：

- `beta_i`
- `beta_total`
- numerator/denominator；
- 统计不确定度；
- 方法定义、能群结构和归一化说明。

### 当前代码流程

当前 `beta_eff.h5` 输出：

- `beta_i`
- `beta_total`
- `uncertainty`，但全为 0；
- `diagnostics/numerator`
- `diagnostics/denominator`
- 材料信息；
- `method_comparison/method_e/*`。

### 差异影响

当前输出适合做 deterministic-style 后处理结果记录，但还不具备完整统计报告能力。

## 11. 主要差异汇总表

| 项目 | 理论计算流程 | 当前代码计算流程 | 影响 |
|---|---|---|---|
| 伴随量 | 真伴随通量或等价重要性 | response-weighted importance | 不能直接套所有严格伴随公式 |
| 能量变量 | 区分入射能量 `E` 和出生能量 `E'` | 主要使用 source group 与 response collision group | 出生谱 `chi` 不宜直接混用 |
| delayed 信息 | 在 `F_d,k` 或 `chi_d,k` 中进入 | 在 transfer function family 分份中进入 | 后处理分子必须乘 `F_total` |
| 分子 | `phi_dagger * F_d,k * phi` | `I_delayed_k * F_total` | delayed yield 已提前计入 |
| 分母 | `phi_dagger * F_t * phi` | `(I_prompt + sum I_delayed) * F_total` | family 求和作为总 importance |
| 能群处理 | 可保留 `I(r,g) * F(r,g)` | 当前 Method E 用 `I(cell)` 乘 `sum_g F_g` | 丢失 family importance 的能群相关性 |
| 空间处理 | 连续空间积分 | mesh cell 求和 + 材料采样均匀化 | 异质区存在离散误差 |
| 统计误差 | ratio estimator | uncertainty 写 0 | 暂无真实不确定度 |
| 源重要性 | 理想 `I*(cell,E)` | `q(cell) * chi_empirical(g|cell)` | 源能群重要性仍是近似 |

## 12. 对后续改造的直接启示

最值得优先处理的不是归一化常数，而是以下三点：

1. **使用 family 分群数据**

当前 `adjoint_flux.h5` 已经写出：

```text
family_resolved/*/flux_group_mean
```

但 Method E 没有使用。应优先改成：

```text
D = sum_cell volume * sum_g I_total(cell,g) * F_total(cell,g)
N_k = sum_cell volume * sum_g I_delayed_k(cell,g) * F_total(cell,g)
```

2. **在文档和输出中明确量的语义**

当前 `adjoint_flux` 不是 true adjoint flux，而是 response-weighted importance。后续公式命名应避免混用 `phi_dagger` 与严格输运伴随概念。

3. **增加 batch-wise 统计**

应保存每个 batch 的：

```text
D_b
N_{k,b}
beta_{k,b}
```

再计算 ratio estimator 的均值和不确定度。

## 13. 结论

理论流程与当前代码流程的最大区别是：理论上 `beta_eff` 是用真正伴随权重直接加权 delayed source 和 total source；当前代码则先在上游 transfer function 中按 prompt/delayed family 拆分贡献，卷积得到 family-resolved importance，再用该 importance 乘 `F_total` 得到分子和分母。

因此当前代码的 Method E 不是严格公式的逐项直接实现，而是一套“delayed 分份前移到 transfer function 阶段”的等价意图实现。它的关键近似和误差来源在于：

- 伴随量不是严格输运伴随通量；
- Method E 后处理把 family importance 折叠成 cell 标量；
- 源能群重要性由经验源谱展开；
- 空间和材料核数据做了 mesh-level 均匀化；
- 尚未实现统计不确定度。

