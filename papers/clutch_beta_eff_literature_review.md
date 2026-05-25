# CLUTCH 与 beta_eff 相关理论调研报告

日期：2026-05-25

## 1. 调研范围与核心结论

本报告围绕项目当前的 `clutch_on`、Green Function、fission matrix、adjoint flux 与 `beta_effective.cpp` 计算链路，重点阅读和对照了以下几类文献：

- β_eff 的经典扰动理论与 Monte Carlo 估计方法；
- IFP/NFP 作为伴随重要性估计的理论基础；
- CLUTCH 方法及其在连续能量 Monte Carlo 灵敏度计算中的用途；
- 裂变矩阵方法生成伴随裂变源的重要性函数；
- OpenMC 中 adjoint-weighted kinetics parameter 与 sensitivity capability 的实现思路。

核心结论如下：

1. `beta_eff` 的严格定义本质上是“缓发裂变中子源对堆系统后续裂变链贡献的重要性”与“总裂变中子源重要性”的比值，而不是简单的核数据比值 `nu_delayed / nu_total`。
2. 伴随通量或伴随重要性函数是 β_eff 的关键权重。文献中常用 IFP/NFP 或 CLUTCH 来避免显式求解连续能量伴随输运方程。
3. IFP 理论上更接近严格伴随重要性，但需要保存祖先/后代关系，内存代价高；NFP 是一代近似，便宜但有截断误差；CLUTCH 用贡献理论和 tracklength/fission-event importance 来降低内存压力，特别适合连续能量 Monte Carlo 的灵敏度计算。
4. 当前项目的 Method E 与 CLUTCH 思想最接近：在上游 Green Function tally 阶段按 `nu_prompt / nu_total`、`nu_delayed,k / nu_total` 分 family 累积重要性，再在 `beta_effective.cpp` 中用 family-resolved importance 乘总裂变源强 `F_total`。这避免了把 `nu_delayed,k` 在分子里重复计数。
5. 当前实现需要特别警惕一个语义边界：`adjoint_flux.h5` 中写明的量是 response-weighted importance field，并不是真正的输运伴随通量 `phi_dagger(r,E)`。因此用它直接做出生谱 `chi_p`、`chi_d,k` 加权会产生物理错位；把 family 选择性放到 `nu_d,k / nu_t` 上更符合当前数据的能量语义。

## 2. beta_eff 的理论定义

反应堆点堆动力学中的有效缓发中子份额可写为各先驱核群贡献之和：

```math
\beta_{\mathrm{eff}} = \sum_i \beta_{i,\mathrm{eff}}
```

严格扰动理论下，若有真正的多群正向通量 `phi_g(r)` 与伴随通量 `phi^\dagger_g(r)`，第 `i` 个缓发群可表示为：

```math
\beta_{i,\mathrm{eff}}
=
\frac{
  \int_V A_{d,i}(r) F_{d,i}(r)\,dr
}{
  \int_V A_p(r) F_t(r)\,dr
}
```

其中：

```math
A_p(r)=\sum_{g'}\phi^\dagger_{g'}(r)\chi_{p,g'}
```

```math
A_{d,i}(r)=\sum_{g'}\phi^\dagger_{g'}(r)\chi_{d,i,g'}
```

```math
F_t(r)=\sum_g \nu_{t,g}\Sigma_{f,g}(r)\phi_g(r)
```

```math
F_{d,i}(r)=\sum_g \nu_{d,i,g}\Sigma_{f,g}(r)\phi_g(r)
```

物理含义是：缓发中子不是只按“产生数量”计入，而是按它们在相空间中引发未来裂变链的能力计入。因此 β_eff 通常不同于简单材料量：

```math
\beta_0 = \frac{\int \nu_d\Sigma_f\phi\,drdE}{\int \nu_t\Sigma_f\phi\,drdE}
```

Meulekamp 与 van der Marck 的 Monte Carlo β_eff 论文强调：经典计算需要伴随和能谱加权，直接显式求伴随函数不方便；通过伴随函数的物理解释，可以在常规 `k_eff` Monte Carlo 计算中构造 β_eff 估计量。

## 3. IFP、NFP 与伴随重要性

### 3.1 IFP 的物理解释

IFP，即 Iterated Fission Probability，把某相空间位置引入一个中子后的“未来裂变链贡献”作为该中子的伴随重要性。Kiedrowski 等的 adjoint-weighted tally 文献以及 OpenMC kinetics 文档都采用这一解释：在前向功率迭代中记录若干代祖先信息，用后续裂变贡献给早先的源/碰撞位置赋权。

OpenMC 官方用户文档中的 IFP 动力学参数实现给出：

```math
\beta_{\mathrm{eff}}
=
\frac{S_{\mathrm{ifp-beta-numerator}}}
       {S_{\mathrm{ifp-denominator}}}
```

`Lambda_eff` 还需要除以 `k_eff`：

```math
\Lambda_{\mathrm{eff}}
=
\frac{S_{\mathrm{ifp-time-numerator}}}
       {S_{\mathrm{ifp-denominator}} k_{\mathrm{eff}}}
```

这说明在 IFP 框架里，β_eff 的核心仍是“带伴随重要性权重的 delayed numerator / denominator”。

### 3.2 NFP 的一代近似

NFP，即 Next Fission Probability，只统计下一代裂变贡献。Meulekamp 与 van der Marck 指出，NFP 不再严格正比于伴随函数，但对于 β_eff 这类全局积分量通常仍有较好精度，因为有效性修正只作用在 `beta_eff` 与基本 delayed fraction 的差异上。代价是存在一代截断误差。

### 3.3 IFP 的工程问题

IFP 需要保存粒子祖先、潜伏代数或 tally 信息。文献中反复指出两个问题：

- 内存压力大，尤其是空间-能量分辨 tally；
- 潜伏代数越大，截断误差越小，但统计噪声和计算代价会上升。

这也是 CLUTCH 和 fission-matrix-based importance 被提出的直接动机。

## 4. CLUTCH 方法理论要点

CLUTCH 全称通常写作 Contribution-Linked eigenvalue sensitivity/Uncertainty estimation via Tracklength importance CHaracterization。它的核心不是显式解伴随输运方程，而是用贡献理论中的“重要性/贡献”解释，在前向 Monte Carlo 过程中构造用于灵敏度和不确定度分析的伴随权重。

从 RMC 的 CLUTCH 文献和 OpenMC sensitivity capability 文献可归纳出以下要点：

1. CLUTCH 目标最初主要是 `k_eff` 对核数据的灵敏度系数，而非专门的 β_eff；但它处理的是同一类伴随加权问题。
2. IFP 理论更直接，但内存压力大；CLUTCH 通过 tracklength/collision/fission-event 贡献量来减少保存历史链的负担。
3. RMC 文献区分了 C-CLUTCH 与 F-CLUTCH：C-CLUTCH 需要保存每次碰撞的相关反应率，F-CLUTCH 只保存裂变点相关反应率，后者内存更低、效率更高。
4. 裂变矩阵方法可以生成 CLUTCH 所需的伴随裂变源分布，即用左特征向量或高阶模式信息近似源重要性。
5. OpenMC sensitivity 文献提出结合 IFP 与 CLUTCH 的统一伴随/广义伴随框架，并指出裂变矩阵方法更适合生成 CLUTCH 的伴随源分布，以规避 IFP 的潜伏代数和统计波动问题。

因此，CLUTCH 对本项目的启发是：如果目标是 β_eff，而不是核数据灵敏度，也可以借鉴它“在前向过程里累积响应贡献、用裂变矩阵生成源端重要性、再做响应侧卷积”的结构。

## 5. 当前项目实现与文献理论的对应

### 5.1 当前主要链路

当前代码大体形成如下链路：

1. `FissionMatrix` 在源状态空间上构造裂变矩阵并求伴随源重要性 `I*(source_state)`。
2. `GreenFunctionMesh` 记录从源状态到响应位置的 transfer function。
3. `AdjointFlux` 做卷积：

```math
\Phi^\dagger(r)
=
\sum_{s} T(s \rightarrow r) I^*(s)
```

1. `tally_scoring.cpp` 在 `SCORE_GREENFUNCTION` 中，按 family 将裂变贡献拆为：

```math
w_{\mathrm{prompt}} = \nu_t \frac{\nu_p}{\nu_t}
```

```math
w_{\mathrm{delayed},k} = \nu_t \frac{\nu_{d,k}}{\nu_t}
```

代码变量里分母是 `nu_total_yield`，物理上就是总产额。

1. `beta_effective.cpp` 的 Method E 用上游 family-resolved importance：

```math
D = \sum_c \Delta V \, I_{\mathrm{total}}(c) F_t(c)
```

```math
N_k = \sum_c \Delta V \, I_{\mathrm{delayed},k}(c) F_t(c)
```

```math
\beta_{k,\mathrm{eff}} = N_k / D
```

这里分子乘 `F_total` 而不是 `F_delayed,k` 是正确的，因为 `I_delayed,k` 在上游 tally 阶段已经包含了 `nu_delayed,k / nu_total` 的分份权重；再乘 `F_delayed,k` 会导致 `nu_delayed,k` 双重计数。

### 5.2 与严格伴随公式的差别

严格公式要求的 `phi^\dagger(r,E)` 是输运伴随通量，其能群索引对应出生谱或输运相空间权重。而当前 `adjoint_flux.cpp` 的元数据明确说明：

- `physical_quantity` 是 response-weighted importance field，不是真正的 adjoint flux；
- group index 表示响应位置导致裂变的碰撞能量，不是裂变中子出生能量。

这会直接影响不同算法公式的可信度：

- 若用当前响应碰撞能量分群的 `Phi^\dagger_g(r)` 去乘 `chi_p,g` 或 `chi_d,k,g`，就是把“碰撞能量”误作“出生能量”，物理语义不一致。
- Method E 把 family 选择性放在 `nu_d,k(E_collision) / nu_total(E_collision)` 上，与当前 group 语义一致，因此比直接 `chi` 加权更稳健。

### 5.3 当前 `clutch.cpp` 的状态

`src/clutch.cpp` 里的 `CLUTCH_TEST()` 当前为空函数，真正承载 CLUTCH-like 数据流的是：

- `SCORE_GREENFUNCTION` 的 family-resolved transfer function tally；
- `FissionMatrix` 的 source-state 伴随源；
- `AdjointFlux` 的 source-state 卷积；
- `BetaEffective` 的 Method E。

因此，项目中的“clutch”更准确地说是一个开关和数据流触发机制，而不是已经独立封装完成的 CLUTCH 算法模块。

## 6. 关键文献阅读摘要

### 6.1 Meulekamp and van der Marck: Calculating the Effective Delayed Neutron Fraction with Monte Carlo

这篇文献直接面向 β_eff。主要贡献是提出可嵌入常规 `k_eff` Monte Carlo 的 β_eff 估计量，不显式求解伴随通量。其理论出发点是 IFP 与伴随函数成比例，工程上又讨论了 NFP 近似。文献指出 NFP 不严格等价于 IFP，但对 β_eff 这类全局量通常可接受。

对本项目的意义：

- β_eff 估计应围绕“未来裂变贡献的重要性”构造；
- 只做材料 delayed fraction 不够；
- 对全局 β_eff，近似伴随函数的误差可能比局部量更容易被积分平均，但逐 delayed group 的分配仍可能敏感。

### 6.2 Kiedrowski, Brown, and Wilson: Adjoint-Weighted Tallies for k-Eigenvalue Calculations with Continuous-Energy Monte Carlo

该文系统说明了在连续能量 `k` 本征值 Monte Carlo 中做 adjoint-weighted tally 的方法。核心做法是用 IFP 估计基本伴随模，并在前向功率迭代中对 tally 贡献加权。文献特别强调潜伏代数、截断误差、统计噪声和内存的折中。

对本项目的意义：

- β_eff、反应性微扰、adjoint-weighted flux 都属于同一类 ratio-of-adjoint-weighted-integrals 问题；
- 增加能量/空间分辨率会迅速推高内存；
- 当前采用 fission matrix + transfer function 的路径，工程上是在绕开完整 IFP 的存储压力。

### 6.3 Calculation of adjoint-weighted reactor kinetics parameters in OpenMC

该文在 OpenMC 中实现了有效缓发中子份额和有效中子代时间，比较了 NFP、IFP 和 CLUTCH 三种伴随权重解释，并用 MCNP5 与实验数据进行基准验证。

对本项目的意义：

- OpenMC 语境下 β_eff 可用多种伴随重要性近似；
- NFP/IFP/CLUTCH 的差别主要在伴随权重如何估计；
- 当前项目若要证明 Method E，应对标该文中的 adjoint-weighted kinetics parameter 框架，而不是只对标材料 β。

### 6.4 Computing eigenvalue sensitivity coefficients based on the CLUTCH method with RMC

该文面向 RMC 中连续能量核数据灵敏度系数。它指出 IFP 理论准确但内存消耗大，提出 C-CLUTCH 与 F-CLUTCH，并用裂变矩阵生成 CLUTCH 所需的伴随裂变源分布。结果显示 F-CLUTCH 与 C-CLUTCH 精度相当，但内存和计算效率更好。

对本项目的意义：

- F-CLUTCH 的“只在裂变点保存贡献”与当前 `SCORE_GREENFUNCTION` 的裂变事件 family 累积更接近；
- 裂变矩阵左特征向量作为伴随源是合理路线；
- 若后续要把 CLUTCH 从 beta_eff 专用流程扩展到灵敏度计算，可优先参考 F-CLUTCH。

### 6.5 Development of continuous-energy sensitivity analysis capability in OpenMC

该文提出 OpenMC 中结合 IFP 与 CLUTCH 的统一伴随/广义伴随框架，用于 `k_eff` 和反应率比对核数据的灵敏度计算。文中指出 IFP 的潜伏代数存在收敛与统计涨落问题，裂变矩阵方法更适合生成 CLUTCH 伴随源分布。

对本项目的意义：

- 当前 source-state fission matrix 改造方向与该文建议一致；
- `I*(cell,g_source)` 比单纯 `I*(cell)` 更适合降低源端能量压缩误差；
- 需要保持 collapsed 路径用于回归验证。

### 6.6 Theory and applications of the fission matrix method for continuous-energy Monte Carlo

该文说明裂变矩阵不仅可估计基本模裂变分布和 dominance ratio，也可提供更高阶前向/伴随裂变源信息。对当前项目而言，它是 `FissionMatrix -> I* -> AdjointFlux` 路线的理论支撑。

## 7. 对当前 beta_eff 方法的判断

当前 Method E 在项目已有方法中理论一致性最好，原因不是数值技巧，而是它与数据语义更匹配：

- 当前 transfer function 的 group/family 信息来自响应侧裂变碰撞；
- 因此 delayed group 选择性应使用 `nu_delayed,k(E_collision) / nu_total(E_collision)`；
- `chi_d,k(E_birth)` 应属于出生谱权重，不能直接拿 response collision energy 的 group 去乘；
- `I_delayed,k * F_total` 避免了 delayed yield 双重计数。

仍然存在的主要误差来源：

1. `I*(source_state)` 当前若由经验 `chi(g|cell)` 展开，仍可能不是严格的能量分辨伴随源。
2. `adjoint_flux` 不是输运伴随通量，不能不加区分地用于所有扰动理论公式。
3. `BetaEffective` 当前读取 family-resolved 数据时，如果只消费 `flux_mean` 而不消费 `flux_group_mean`，会丢失 family 内部能量结构。
4. 网格单元内材料采样、混合材料核数据加权、参考能量下抽样 `chi` 等会引入第二层近似。

## 8. 后续建议

建议按以下顺序推进，而不是直接继续调整 β_eff 归一化常数：

1. 明确术语：在文档和 HDF5 元数据中避免把 response-weighted importance 直接称为 true adjoint flux。
2. 保留 Method E 作为主结果路径，并在报告/输出中说明其 family-resolved 权重已经包含 `nu_d,k / nu_t`。
3. 检查 `BetaEffective` 是否充分使用 `family_resolved/*/flux_group_mean`。如果只用 `flux_mean`，应增加分群 Method E 版本：

```math
N_k =
\sum_c \Delta V
\sum_g I_{\mathrm{delayed},k}(c,g)
      \nu_{t,g}\Sigma_{f,g}(c)\phi_g(c)
```

```math
D =
\sum_c \Delta V
\sum_g I_{\mathrm{total}}(c,g)
      \nu_{t,g}\Sigma_{f,g}(c)\phi_g(c)
```

1. 建立三个验证层级：

- 均匀单材料体系：验证 `sum beta_i` 与核数据 delayed fraction 的量级一致；
- ICSBEP/MCNP/Serpent benchmark：验证总 β_eff；
- delayed group-wise benchmark：验证各群分配，重点观察 Method E 是否修复 Method B 的 family redistribution。

1. 若继续靠近严格扰动理论，应把源端从 `I*(cell)` 稳定升级为 `I*(cell,g_source)`，并确保 Green Function 和 AdjointFlux 都使用同一 source-state 索引。

## 9. 已阅读或纳入的本地文献

- `papers/Calculating the Effective Delayed Neutron Fraction with Monte Carlo.pdf`
- `papers/Calculation of adjoint-weighted reactor kinetics parameters in OpenMC.pdf`
- `papers/Adjoint-Weighted Tallies for k-Eigenvalue Calculations with Continuous-Energy Monte Carlo.pdf`
- `papers/Computing eigenvalue sensitivity coefficients to nuclear data based on the CLUTCH method with RMC code.pdf`
- `papers/Development of continuous-energy sensitivity analysis capability in OpenMC.pdf`
- `papers/Theory and applications of the fission matrix method for continuous-energy Monte Carlo.pdf`
- `papers/Monte Carlo and deterministic computational methods for the calculation of the effective delayed neutron fraction.pdf`
- `docs/source/usersguide/kinetics.rst`
- `markdown/分析0423.md`
- `markdown/能量伴随源分布计算修改记录.md`

以下 PDF 本轮未能可靠抽取文本，后续建议使用 OCR 或专用 PDF 工具重新处理：

- `papers/Improvement of the CLUTCH method for sensitivity analysis of k-eigenvalue to continuous-energy nuclear data in NECP-MCX.pdf`
- `papers/Proposal_of_direct_calculation_of_kinetic_parameters_beta_sub_eff_and_Lambda_based_on_continuous_ene1.pdf`
