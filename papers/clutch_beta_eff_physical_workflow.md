# 基于 CLUTCH 思想的 beta_eff 物理计算流程

日期：2026-05-25

本文只从文献方法和物理定义出发整理，不引用当前项目代码实现。

## 1. 物理目标

有效缓发中子份额 `beta_eff` 衡量的是：缓发中子源对维持临界裂变链的有效贡献，占全部裂变中子源有效贡献的比例。

它不是单纯的核数据比例：

```math
\beta_0 = \frac{\nu_d}{\nu_t}
```

而是伴随重要性加权后的比例：

```math
\beta_{\mathrm{eff}}
=
\frac{\langle \phi^\dagger, F_d \phi \rangle}
     {\langle \phi^\dagger, F_t \phi \rangle}
```

其中：

- `phi` 是正向中子通量；
- `phi^\dagger` 是伴随重要性函数；
- `F_d` 是缓发裂变中子源算符；
- `F_t = F_p + F_d` 是总裂变中子源算符；
- 尖括号表示对空间、能量、方向的积分。

分群形式为：

```math
\beta_{\mathrm{eff}}=\sum_k \beta_{k,\mathrm{eff}}
```

```math
\beta_{k,\mathrm{eff}}
=
\frac{\langle \phi^\dagger, F_{d,k}\phi \rangle}
     {\langle \phi^\dagger, F_t\phi \rangle}
```

## 2. CLUTCH 在这里提供什么

CLUTCH 的核心作用是提供一个可在前向 Monte Carlo 计算中估计的伴随重要性权重。

从物理上说，一个中子在某个相空间状态 `x = (r,E,Omega)` 的伴随重要性，可理解为它对未来裂变链或指定响应的贡献能力。IFP 用未来若干代裂变后代数来估计这个重要性；CLUTCH 则用贡献理论，把当前粒子轨迹、碰撞或裂变事件对响应的贡献与源重要性联系起来，避免完整保存长链祖先信息。

因此，基于 CLUTCH 计算 `beta_eff` 的关键不是先显式求出完整 `phi^\dagger`，而是在 Monte Carlo 过程中构造等价的伴随加权积分：

```math
N_k \approx \text{CLUTCH-weighted delayed group k fission source}
```

```math
D \approx \text{CLUTCH-weighted total fission source}
```

最后：

```math
\beta_{k,\mathrm{eff}} = N_k / D
```

## 3. 连续能量物理表达式

对第 `k` 个缓发群，严格表达式可写为：

```math
N_k =
\int d\mathbf r
\int dE
\int dE'
\phi(\mathbf r,E)
\Sigma_f(\mathbf r,E)
\nu_{d,k}(\mathbf r,E)
\chi_{d,k}(E')
\phi^\dagger(\mathbf r,E')
```

分母为：

```math
D =
\int d\mathbf r
\int dE
\int dE'
\phi(\mathbf r,E)
\Sigma_f(\mathbf r,E)
\nu_t(\mathbf r,E)
\chi_t(E')
\phi^\dagger(\mathbf r,E')
```

其中：

- `E` 是诱发裂变的入射中子能量；
- `E'` 是裂变产生中子的出生能量；
- `nu_{d,k}` 是第 `k` 个缓发群平均产额；
- `chi_{d,k}` 是第 `k` 个缓发群出生谱；
- `nu_t` 与 `chi_t` 对应总裂变中子源。

如果采用 prompt/delayed 分解，则：

```math
F_t = F_p + \sum_k F_{d,k}
```

## 4. 基于 CLUTCH 的计算流程

### 步骤 1：完成临界正向 Monte Carlo 计算

先进行常规 `k_eff` 本征值计算，得到稳定的正向裂变源分布和正向通量分布。

这一阶段提供：

- 空间-能量上的裂变反应率；
- 每次裂变事件的入射能量、位置、材料或核素；
- 裂变中子总产额 `nu_t(E)`；
- 缓发群产额 `nu_{d,k}(E)`；
- 必要时的出生谱 `chi_p(E')`、`chi_{d,k}(E')`。

### 步骤 2：定义 CLUTCH 响应

对 `beta_eff` 来说，需要两个响应：

```math
R_D = \langle \phi^\dagger, F_t\phi \rangle
```

```math
R_{N,k} = \langle \phi^\dagger, F_{d,k}\phi \rangle
```

`R_D` 是总裂变中子源的有效贡献，`R_{N,k}` 是第 `k` 个缓发群裂变中子源的有效贡献。

### 步骤 3：估计源端伴随重要性

CLUTCH 需要某种源端重要性 `I^*(s)`，其中 `s` 是源状态，例如：

```math
s = (\mathbf r_s, E_s)
```

文献中常见做法包括：

- IFP 估计：用未来若干代裂变后代数估计 `I^*`；
- 裂变矩阵估计：构造裂变源状态之间的转移矩阵，取左本征向量作为伴随裂变源重要性；
- CLUTCH/贡献理论估计：用轨迹或裂变事件贡献与源状态建立联系。

物理要求是：`I^*(s)` 应表示一个出生于源状态 `s` 的中子对未来临界裂变链的相对贡献。

### 步骤 4：沿粒子历史建立贡献链接

对每个源中子及其后代历史，记录该历史在不同位置和能量发生裂变事件时，对响应的贡献。

对于一次由入射能量 `E` 在位置 `r` 诱发的裂变，总裂变贡献可写为：

```math
w_t(r,E) \propto \nu_t(E)\Sigma_f(r,E)
```

第 `k` 个缓发群贡献为：

```math
w_{d,k}(r,E) \propto \nu_{d,k}(E)\Sigma_f(r,E)
```

如果在贡献 tally 阶段先按产额分份拆分，则可写为：

```math
w_{d,k}(r,E)
=
w_t(r,E)\frac{\nu_{d,k}(E)}{\nu_t(E)}
```

这一步就是把 delayed family 的选择性放在诱发裂变的入射能量侧。

### 步骤 5：用 CLUTCH 重要性加权贡献

每个贡献事件需要乘以其所属源状态的重要性 `I^*(s)`。

总分母累积：

```math
D
\leftarrow
D + I^*(s)\, w_t(r,E)
```

第 `k` 群分子累积：

```math
N_k
\leftarrow
N_k + I^*(s)\, w_{d,k}(r,E)
```

若保留出生谱加权，则 delayed 分子还需要包含 `chi_{d,k}` 与伴随重要性在出生能量上的折叠；若采用 CLUTCH 的 family 分份形式，则 `nu_{d,k}/nu_t` 已经承担 delayed group 选择性，后处理时不能再次乘 `nu_{d,k}`。

### 步骤 6：批次平均和归一化

对所有 active batches 统计：

```math
\bar D = \frac{1}{B}\sum_b D_b
```

```math
\bar N_k = \frac{1}{B}\sum_b N_{k,b}
```

然后计算：

```math
\beta_{k,\mathrm{eff}}
=
\frac{\bar N_k}{\bar D}
```

```math
\beta_{\mathrm{eff}}
=
\sum_k \beta_{k,\mathrm{eff}}
```

统计不确定度应按 ratio estimator 处理，而不是简单分别给 `N_k` 和 `D` 求相对误差。

## 5. 推荐的分群计算形式

如果用多群离散形式，可写为：

```math
D =
\sum_c \Delta V_c
\sum_g
I_t^*(c,g)\,
\nu_{t,g}(c)\Sigma_{f,g}(c)\phi_g(c)
```

第 `k` 群：

```math
N_k =
\sum_c \Delta V_c
\sum_g
I_{d,k}^*(c,g)\,
\nu_{t,g}(c)\Sigma_{f,g}(c)\phi_g(c)
```

这里 `I_{d,k}^*` 已经是按 delayed group `k` 分份后的 CLUTCH 重要性。如果没有提前分份，则应使用：

```math
N_k =
\sum_c \Delta V_c
\sum_g
I^*(c,g)\,
\nu_{d,k,g}(c)\Sigma_{f,g}(c)\phi_g(c)
```

两种写法不能混用：

- 若 `I_{d,k}^*` 已含 `nu_{d,k}/nu_t`，分子乘 `nu_t Sigma_f phi`；
- 若 `I^*` 未含 delayed group 信息，分子乘 `nu_{d,k} Sigma_f phi`；
- 不能同时用 `I_{d,k}^*` 和 `nu_{d,k}`，否则 delayed yield 双重计数。

## 6. CLUTCH-beta_eff 与 IFP-beta_eff 的关系

IFP-beta_eff 可理解为直接用未来裂变后代数估计 `phi^\dagger`，再计算 delayed source 与 total source 的伴随加权比值。

CLUTCH-beta_eff 则是用贡献链接方式重写同一个伴随加权比值：

```math
\frac{\langle \phi^\dagger, F_{d,k}\phi \rangle}
     {\langle \phi^\dagger, F_t\phi \rangle}
\quad
\Longrightarrow
\quad
\frac{\text{source-importance-weighted delayed contribution}}
     {\text{source-importance-weighted total contribution}}
```

因此两者目标相同，差别在于伴随重要性和加权积分的估计方式：

- IFP 更直接、更接近定义，但内存和潜伏代数代价高；
- CLUTCH 更工程化，用贡献链接减少长历史存储，适合连续能量大规模计算；
- 裂变矩阵可作为 CLUTCH 源端重要性的稳定估计器。

## 7. 物理一致性检查

整理 CLUTCH-beta_eff 流程时，需要检查以下物理一致性：

1. 分子和分母使用同一个伴随重要性定义。
2. delayed group 选择性只引入一次。
3. 如果使用 `chi_d,k`，其能量变量必须是裂变中子出生能量，而不是诱发裂变的入射能量。
4. 如果使用 `nu_{d,k}/nu_t`，其能量变量是诱发裂变的入射能量。
5. `beta_eff` 是 ratio estimator，不能把分子、分母分别归一化到不一致的源强。
6. `sum_k beta_{k,eff}` 应为总 `beta_eff`，量级应接近但不必等于材料 delayed fraction。

## 8. 最简流程图

```text
正向 k_eff Monte Carlo
        |
        v
统计裂变事件与源状态
        |
        v
估计源状态伴随重要性 I*(s)
        |
        v
沿历史建立贡献链接：s -> 裂变事件(r,E)
        |
        v
累积分母 D: I*(s) * total fission contribution
        |
        v
累积分子 N_k: I*(s) * delayed-k contribution
        |
        v
beta_k,eff = N_k / D
        |
        v
beta_eff = sum_k beta_k,eff
```

## 9. 结论

基于 CLUTCH 思想计算 `beta_eff`，本质上是把严格定义中的伴随加权积分改写为前向 Monte Carlo 中可统计的贡献链接估计。流程核心是：

1. 用裂变链贡献估计源状态重要性；
2. 用该重要性加权每次裂变事件对总源和缓发群源的贡献；
3. 对 delayed numerator 和 total denominator 做同一套统计平均；
4. 取比值得到各群 `beta_k,eff` 和总 `beta_eff`。

最重要的物理边界是：delayed group 信息只能引入一次；若在 CLUTCH 贡献阶段已经按 `nu_{d,k}/nu_t` 分份，后续分子就必须乘总裂变源，而不是再次乘 delayed 裂变源。

