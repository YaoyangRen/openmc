# CLUTCH beta_eff 改造软件报告

日期：2026-06-11

## 1. 报告范围

本文档说明当前项目中 `beta_eff` 计算相关改造的实现状态和运行流程。重点覆盖：

- `cell x material` 源状态定义。
- F-CLUTCH 主方法及 missing-I fallback。
- C-CLUTCH transfer-function 方法。
- 新增 IFP-ancestry 诊断方法。
- 新增 CLUTCH-IFP 方法，即使用 IFP 得到的重要性函数替换 fission-matrix `I*`。
- 从一个 k-eigenvalue 算例开始输运，到输出 `beta_eff.h5`、`generation_time.h5`、`fission_matrix.h5` 的完整流程。

本文档描述的是当前代码实现，不引入出生能量维度。

## 2. 改造目标

原始 material-resolved CLUTCH 使用 inactive batches 中统计的 fission matrix 求解源重要性：

```text
I*(cell, material)
```

然后在 active batches 中用 fission site 或 transfer function 折叠 delayed group 贡献。实际大网格、多材料算例中，F-CLUTCH 会遇到一部分 active fission site 的 `I*(cell,material)` 缺失，导致 site 被 drop，结果偏低。

本次改造的目标是：

1. 保留 F-CLUTCH 为根级主输出。
2. 给 F-CLUTCH 增加 missing-I fallback，减少 active fission site 丢失。
3. 新增 IFP-ancestry 诊断方法，用 N 代祖先 delayed group 直接估计 beta_eff。
4. 新增 CLUTCH-IFP 方法，用 IFP ancestry 得到的 `I_ifp(cell,material)` 替代 fission-matrix `I*`，作为新的并列方法输出。
5. 不修改出生能量维度，不引入核素维度。

## 3. 主要代码改动

### 3.1 新增 CLUTCH-IFP ancestry bank

新增文件：

- `include/openmc/clutch_ifp.h`
- `src/clutch_ifp.cpp`

新增全局 bank 数据：

```cpp
simulation::clutch_ifp_source_state_bank
simulation::clutch_ifp_source_delayed_group_bank
simulation::clutch_ifp_fission_state_bank
simulation::clutch_ifp_fission_delayed_group_bank
```

它们的作用与 OpenMC 原生 IFP bank 类似，但额外携带 `cell x material` source state。OpenMC 原生 IFP bank 只携带 delayed group 和 lifetime，不能直接得到 material-resolved source-state importance，因此这里新增 CLUTCH 专用 ancestry bank。

### 3.2 新增 CLUTCH-IFP bank 搬运流程

改动文件：

- `src/bank.cpp`
- `src/eigenvalue.cpp`
- `src/simulation.cpp`

实现内容：

- fission bank 排序时同步排序 CLUTCH-IFP ancestry bank。
- fission bank 抽样成下一代 source bank 时，同步复制 CLUTCH-IFP ancestry 数据。
- MPI 分发 source bank 时，同步发送、接收、反序列化 CLUTCH-IFP ancestry 数据。
- 串行路径下直接复制 temporary ancestry bank 到 source ancestry bank。

### 3.3 CE/MG 裂变点接入 ancestry 记录

改动文件：

- `src/physics.cpp`
- `src/physics_mg.cpp`

每个成功写入 fission bank 的 fission site 都会追加：

```text
child_state = cell(site.r) x material(current collision material)
child_delayed_group = site.delayed_group
```

这样下一代 source neutron 能携带最近 N 代祖先的 source state 和 delayed group。

### 3.4 BetaEffectiveAccumulator 扩展

改动文件：

- `include/openmc/beta_effective_accumulator.h`
- `src/beta_effective_accumulator.cpp`

新增内容：

- `score_ifp_ancestry_event(...)`
- `compute_ifp_ancestry_result()`
- `compute_clutch_ifp_result()`
- `fold_current_clutch_ifp_batch()`
- `write_clutch_ifp_method_group(...)`
- F-CLUTCH missing-I fallback。
- `beta_eff.h5` 中新增 `method/clutch_ifp` 和相关 diagnostics。

## 4. 源状态定义

当前源状态固定为：

```text
state = cell * n_materials + material_index
```

其中：

- `cell` 来自动力学 mesh，即 `kinetics_mesh`。
- `material_index` 使用 OpenMC 内部材料下标。
- `n_materials = model::materials.size()`。

当前方法不使用：

- 出生能量。
- 入射能量。
- delayed family 作为 source-state 维度。
- 核素维度。
- 角度维度。

## 5. 从算例输运开始的计算流程

本节以一个普通 eigenvalue 算例为例，例如 ZPPR20C 类似输入：

```xml
<run_mode>eigenvalue</run_mode>
<particles>...</particles>
<batches>...</batches>
<inactive>...</inactive>
<kinetics_mesh>...</kinetics_mesh>
<adjoint_source>...</adjoint_source>
```

### 5.1 初始化阶段

OpenMC 初始化几何、材料、截面、settings 后，当前项目额外初始化以下对象：

1. Shared kinetics mesh

```text
SharedMeshGrid
```

这个网格同时用于：

- fission matrix source state。
- beta_eff source-state scoring。
- generation time denominator。
- sensitivity 中 collapsed cell importance。

2. FissionMatrix

当 `settings::clutch_on` 开启时创建：

```text
FissionMatrix(grid, batches, kinetics_energy_edges, score_start_batch)
```

虽然构造函数仍保留 `kinetics_energy_edges` 参数，但当前 `beta_eff` 主流程不使用出生能量维度。

3. BetaEffectiveAccumulator

当 `settings::beta_effective_on` 开启时创建：

```text
BetaEffectiveAccumulator(grid)
```

4. CLUTCH-IFP ancestry banks

当满足：

```text
beta_effective_on
run_mode == eigenvalue
solver_type == MONTE_CARLO
```

时分配：

```text
clutch_ifp_source_state_bank
clutch_ifp_source_delayed_group_bank
clutch_ifp_fission_state_bank
clutch_ifp_fission_delayed_group_bank
```

其 ancestry 长度使用：

```text
N_gen = settings::ifp_n_generation, if set
N_gen = DEFAULT_IFP_N_GENERATION, otherwise
N_gen <= n_inactive
```

### 5.2 每个 batch 开始

进入 `initialize_batch()` 后：

1. 若当前 batch 是 inactive batch，则 fission matrix 开始新 batch。
2. 若当前 batch 是 active batch，则 `BetaEffectiveAccumulator::begin_batch(batch_id)` 开始 active beta_eff 计分。

active batch 内会清空当前 batch 的临时容器：

```text
current_source_states_
current_source_counts_
current_clutch_ifp_source_counts_
current_fission_site_total_by_state_
current_fission_site_delayed_by_state_
current_clutch_ifp_response_by_state_
current_cclutch_transfer_*
```

### 5.3 初始化每条源粒子历史

每条 source neutron 从 source bank 取出后，OpenMC 执行几何定位和材料定位。

当前项目额外执行：

1. inactive batch 中记录 fission matrix 的 parent source state：

```text
FissionMatrix::record_source_birth(r, source_particle_id, material, weight)
```

2. active batch 中记录 beta_eff active source state：

```text
BetaEffectiveAccumulator::record_source_birth(
    r,
    source_particle_id,
    material,
    source_bank_index
)
```

这里 `source_bank_index = p.current_work() - 1`。它用于访问当前 source neutron 携带的 N 代 ancestry bank。

`record_source_birth()` 同时做两件事：

- 给 C-CLUTCH 记录当前 active source particle 的 `source_particle_id -> state`。
- 给 CLUTCH-IFP 记录 N 代祖先 state 的 source count：

```text
C_s = count(active source particles whose N-gen ancestor state is s)
```

### 5.4 粒子输运与裂变事件

粒子输运过程仍由 OpenMC 原有 CE/MG physics 完成。发生裂变时，当前项目在两个层面计分。

#### 5.4.1 C-CLUTCH collision estimator

在 active batch，裂变碰撞处先计算期望裂变源贡献：

```text
T_total
T_delayed,k
```

CE 路径中使用 nuclide delayed yield：

```text
T_delayed,k = nu_t * nu_delayed,k / nu_total
```

MG 路径中使用 MG delayed-nu-fission：

```text
T_delayed,k = nu_t * delayed_nu_fission_k / nu_fission
```

随后调用：

```text
score_cclutch_fission_event(...)
```

C-CLUTCH 先按当前 source state 暂存 transfer response，在 batch 结束时再乘以 `I*` 折叠。

#### 5.4.2 抽样 fission source site

OpenMC 按 `nu_t` 抽样产生 fission source site。每个 site 包含：

```text
site.r
site.wgt
site.E
site.delayed_group
site.parent_id
site.progeny_id
```

其中：

```text
site.delayed_group = 0     prompt
site.delayed_group = 1..8  delayed group
```

CE 路径由 `sample_fission_neutron()` 抽样 delayed group 和出射能量。

MG 路径由 `sample_fission_energy()` 抽样 outgoing group 和 delayed group。

#### 5.4.3 记录 CLUTCH-IFP fission ancestry

当 fission site 成功写入 fission bank 后，调用：

```text
record_clutch_ifp_fission_site(
    p,
    child_state,
    site.delayed_group,
    fission_bank_index
)
```

其中：

```text
child_state = cell(site.r) * n_materials + p.material()
```

该函数将当前 fission site 的 `child_state` 和 `delayed_group` 追加到 ancestry 链。若 ancestry 链长度超过 `N_gen`，则移除最老项，保留最近 N 代。

#### 5.4.4 F-CLUTCH fission-site scoring

active batch 中每个 fission site 调用：

```text
score_fission_site(site.r, site.wgt, material, delayed_group, lifetime, delay)
```

它记录三类信息：

1. F-CLUTCH 主方法的 I* 加权计分。
2. raw fission-site totals by state，供 CLUTCH-IFP 后续折叠。
3. delayed group raw totals by state，供 CLUTCH-IFP 后续折叠。

raw site totals 为：

```text
F_s     = sum_sites_in_state_s w_site
F_k,s   = sum_delayed_group_k_sites_in_state_s w_site
```

F-CLUTCH 主方法直接使用 fission-matrix `I*(state)`：

```text
D_F,b   = sum_sites w_site * I*(state_site)
N_F,k,b = sum_delayed_group_k_sites w_site * I*(state_site)
```

### 5.5 F-CLUTCH missing-I fallback

若 active fission site 的完整 `I*(cell,material)` 不存在，程序不再立即 drop。当前 fallback 为：

```text
I_fallback(cell) =
    mean over available material states in same cell of I*(cell, material)
```

即：

```text
I_used(cell, material) =
    I*(cell, material), if available
    I_fallback(cell), otherwise
```

只有当同一个 cell 没有任何 material-state importance 时，才 drop 该 fission site。

新增 diagnostics：

```text
total_fallback_sites
total_dropped_sites
```

解释：

- `total_fallback_sites` 表示 missing material-state `I*` 但通过同 cell fallback 成功计分的 site。
- `total_dropped_sites` 表示 fallback 也失败后真正被丢弃的 site。

### 5.6 batch 结束时的折叠

每个 active batch 结束时调用：

```text
fold_current_cclutch_batch()
fold_current_clutch_ifp_batch()
```

#### 5.6.1 C-CLUTCH 折叠

C-CLUTCH 使用当前 batch 中每个 source state 的 transfer response：

```text
D_C,b   = sum_s I*(s) / source_count(s) * T_total(s)
N_C,k,b = sum_s I*(s) / source_count(s) * T_delayed,k(s)
```

其中 `source_count(s)` 是当前 batch 中 active source birth 的 state 计数。

#### 5.6.2 CLUTCH-IFP 重要性构造

CLUTCH-IFP 使用 active batch 内的 N 代 ancestry 信息构造：

```text
R_s = sum_current_fission_events p.wgt_last
      for events whose N-gen ancestor state is s

C_s = count_current_source_particles
      whose N-gen ancestor state is s

I_ifp(s) = R_s / C_s
```

`R_s` 来自：

```text
score_ifp_ancestry_event(p.wgt_last, source_bank_index)
```

`C_s` 来自：

```text
record_source_birth(..., source_bank_index)
```

然后 CLUTCH-IFP 把 active fission-site raw totals by state 与 `I_ifp` 折叠：

```text
D_CI,b   = sum_s F_s   * I_ifp(s)
N_CI,k,b = sum_s F_k,s * I_ifp(s)
```

最终：

```text
beta_CI,k = mean_b(N_CI,k,b) / mean_b(D_CI,b)
```

### 5.7 IFP-ancestry 直接诊断

IFP-ancestry 是直接 N-generation estimator，不使用 `I*`，也不折叠 raw fission sites：

```text
D_IFP,b   = sum_current_fission_events p.wgt_last
N_IFP,k,b = sum_current_fission_events p.wgt_last
            * 1[ancestor_delayed_group == k]
```

最终：

```text
beta_IFP,k = mean_b(N_IFP,k,b) / mean_b(D_IFP,b)
```

该方法用于诊断 ancestry 链和 delayed group 标记是否正常。由于 delayed source 是 analog 抽样，单个算例中它的统计误差可能明显大于 C-CLUTCH。

### 5.8 fission bank 到 source bank 的跨代传播

每一代结束后，OpenMC 对 fission bank 做：

1. 排序。
2. 抽样。
3. MPI 或串行 source bank 同步。

当前项目同步执行 CLUTCH-IFP ancestry bank 的相同操作：

```text
fission ancestry bank -> temp ancestry bank -> source ancestry bank
```

因此每个 source neutron 在进入下一代输运时，都携带一条长度最多为 `N_gen` 的 ancestry 链：

```text
states:         [state_{oldest}, ..., state_{newest}]
delayed_groups: [dg_{oldest},    ..., dg_{newest}]
```

当链长度达到 `N_gen` 时，`front()` 就是 N 代祖先信息。

## 6. inactive batches 中的 fission matrix 与 I*

inactive batches 用于统计 fission matrix：

```text
F(parent_state, child_state)
```

每个 source birth 记录 parent state，每个子 fission site 记录 child state。

inactive 结束时：

1. 合并最后一个 inactive batch。
2. 对每个 parent state 用 source count 归一化：

```text
F_norm(parent, child) = fission_weight(parent, child) / source_count(parent)
```

3. 求解伴随源重要性：

```text
I*(parent) = (1/k) * sum_child F_norm(parent, child) * I*(child)
```

4. 得到 sparse `I*(cell,material)`。
5. 将 `I*` 传给 `BetaEffectiveAccumulator`。
6. 写出 `fission_matrix.h5`。

## 7. active batches 中的 beta_eff 方法汇总

当前 `beta_eff.h5` 中包含四类方法。

### 7.1 F-CLUTCH, 根级主方法

根级输出仍来自 F-CLUTCH：

```text
/beta_i
/beta_total
/uncertainty
```

公式：

```text
D_F   = mean_b sum_sites w_site * I_used(state_site)
N_F,k = mean_b sum_delayed_k_sites w_site * I_used(state_site)
beta_F,k = N_F,k / D_F
```

其中 `I_used` 包含 missing-I fallback。

### 7.2 C-CLUTCH

输出路径：

```text
/method/cclutch
```

公式：

```text
D_C   = mean_b sum_s I*(s)/source_count(s) * T_total(s)
N_C,k = mean_b sum_s I*(s)/source_count(s) * T_delayed,k(s)
beta_C,k = N_C,k / D_C
```

C-CLUTCH 使用期望 delayed production，通常方差较低。

### 7.3 CLUTCH-IFP

输出路径：

```text
/method/clutch_ifp
```

公式：

```text
I_ifp(s) = R_s / C_s

D_CI   = mean_b sum_s F_s   * I_ifp(s)
N_CI,k = mean_b sum_s F_k,s * I_ifp(s)
beta_CI,k = N_CI,k / D_CI
```

它的目标是判断 fission-matrix `I*` 是否是误差来源。

### 7.4 IFP-ancestry

输出路径：

```text
/method/ifp_ancestry
```

公式：

```text
D_IFP   = mean_b sum_events w_event
N_IFP,k = mean_b sum_events w_event * 1[ancestor_dg == k]
beta_IFP,k = N_IFP,k / D_IFP
```

它是直接诊断 estimator，不替代主输出。

## 8. 输出文件

### 8.1 fission_matrix.h5

主要内容：

- sparse COO fission matrix。
- source state metadata。
- material IDs。
- sparse adjoint source `I*(cell,material)`。
- coverage diagnostics。

关键 metadata：

```text
source_state_definition = cell_material
state_indexing = state=cell*n_materials+material_index
matrix_semantics = row=(parent_cell,parent_material), col=(child_cell,child_material)
```

### 8.2 beta_eff.h5

当前版本：

```text
version = 4.1
```

根级：

```text
/beta_i
/beta_total
/uncertainty
```

仍为 F-CLUTCH 主结果。

方法组：

```text
/method/fclutch_cell_material
/method/cclutch
/method/clutch_ifp
/method/ifp_ancestry
```

diagnostics：

```text
total_fission_sites
total_scored_sites
total_dropped_sites
total_fallback_sites
total_invalid_delayed_group_sites
total_cclutch_events
total_cclutch_scored_events
total_ifp_ancestry_events
total_ifp_ancestry_scored_events
total_ifp_ancestry_incomplete_events
```

### 8.3 generation_time.h5

当前 generation time 仍使用 CLUTCH-weighted denominator：

```text
/method/fclutch_cell_material
/method/cclutch
```

本次没有把 CLUTCH-IFP 扩展到 generation time。

## 9. 当前程序终端输出

`beta_eff` 汇总表包含：

```text
group
F-CLUTCH
unc_F
C-CLUTCH
unc_C
CLUTCH-IFP
unc_CI
IFP-ancestry
unc_IFP
```

说明：

- F-CLUTCH 是主方法。
- C-CLUTCH 是低方差 transfer-function 方法。
- CLUTCH-IFP 是用 IFP-derived importance 替换 `I*` 的新增方法。
- IFP-ancestry 是直接 N-generation 诊断方法。

## 10. 与 OpenMC 原生 IFP 的关系

OpenMC 原生 IFP 已有 delayed group ancestry 和 lifetime ancestry，用于 tally：

```text
ifp-beta-numerator
ifp-denominator
ifp-time-numerator
```

但原生 IFP bank 不携带 `cell x material` state，因此无法直接构造：

```text
I_ifp(cell, material)
```

本项目新增的 CLUTCH-IFP bank 复用了 OpenMC IFP 的核心思想：

```text
用 N 代后代裂变概率代表祖先 source state 重要性
```

但扩展了 ancestry 内容：

```text
delayed_group -> delayed_group + source_state(cell,material)
```

因此它可以同时支持：

- IFP-ancestry 直接 beta_eff 诊断。
- CLUTCH-IFP source-state importance。

## 11. 结果分析建议

运行 ZPPR20C 后建议重点检查：

1. F-CLUTCH coverage

```text
total_fission_sites
total_scored_sites
total_fallback_sites
total_dropped_sites
```

若 `total_fallback_sites` 很大，说明 fission-matrix `I*(cell,material)` 覆盖仍不足。

2. F-CLUTCH 与 C-CLUTCH 差异

若 C-CLUTCH 稳定接近参考值，而 F-CLUTCH 仍偏低，说明 fission-site analog 计分和 missing-I coverage 是主要误差来源。

3. CLUTCH-IFP 与 F-CLUTCH 差异

若 CLUTCH-IFP 接近 OpenMC IFP 或 C-CLUTCH，而 F-CLUTCH 偏离，说明 fission-matrix `I*` 可能是主要问题。

4. IFP-ancestry uncertainty

若 `unc_IFP` 较大，说明 direct analog ancestry 的 delayed source 样本仍不足。它更适合作为趋势诊断，不适合作为低方差主结果。

## 12. 验证状态

本次文档对应的源码改动已经完成静态检查：

```text
git diff --check
```

未执行编译和单元测试，因为当前 ZPPR20C 算例正在运行。

当前算例完成后建议执行：

```powershell
cmake --build .\build --target test_beta_effective_accumulator test_fission_matrix openmc
.\build\bin\test_beta_effective_accumulator.exe
.\build\bin\test_fission_matrix.exe
```

## 13. 当前限制

1. CLUTCH-IFP 目前只用于 beta_eff，不用于 generation time。
2. CLUTCH-IFP 的 `I_ifp` 当前按 active batch 构造和折叠，统计质量依赖 active fission event 数量。
3. F-CLUTCH fallback 是同 cell material-state 平均值，不是严格理论解，只是减少 missing-I drop 的工程修正。
4. 未引入出生能量维度，当前判断是 birth-energy source state 对本问题不具备统计收益。
5. OpenMC 原生 IFP tally 仍可作为外部参考，但当前 CLUTCH-IFP 使用的是项目内新增 ancestry bank，因为需要 `cell x material` state。

## 14. 小结

当前程序形成了四条并列 beta_eff 路径：

```text
F-CLUTCH       = 主结果，使用 fission-matrix I*，带 missing-I fallback
C-CLUTCH       = transfer-function 折叠，低方差对照
CLUTCH-IFP     = 用 IFP-derived I_ifp 替换 I* 的新增方法
IFP-ancestry   = 直接 N-generation ancestry 诊断
```

这一结构可以把误差来源拆开：

- 若 F-CLUTCH 与 C-CLUTCH 差异大，先看 missing-I coverage。
- 若 CLUTCH-IFP 改善明显，说明 fission-matrix `I*` 是关键误差来源。
- 若 IFP-ancestry 与 OpenMC 原生 IFP 一致，说明 delayed group ancestry 传播和 delayed 标记逻辑是正确的。
