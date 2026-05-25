# beta_eff 代码计算流程整理

日期：2026-05-25

本文只依据当前代码实现整理，不引用文献理论。涉及的主文件包括：

- `src/simulation.cpp`
- `src/particle.cpp`
- `src/tallies/tally_scoring.cpp`
- `src/fission_matrix.cpp`
- `src/greenfunction_mesh.cpp`
- `src/adjoint_flux.cpp`
- `src/beta_effective.cpp`
- `include/openmc/beta_effective.h`

## 1. 总体执行链路

当前 `beta_eff` 计算不是单独运行的模块，而是在 `clutch_on` 打开后，随本征值计算最后一个 batch 自动触发。

主链路如下：

```text
粒子输运过程中
  |
  |-- FluxMesh 累积正向 track-length 通量
  |
  |-- FissionMatrix 在 inactive batches 统计裂变矩阵
  |
  |-- GreenFunctionMesh 在 fission tally 中统计 source_state -> response cell 传递函数
  v
inactive batches 结束
  |
  |-- FissionMatrix::compute_adjoint_source()
  |-- 写出 fission_matrix.h5
  v
最后一个 batch 结束
  |
  |-- GreenFunctionMesh::finalize_greenfunction_mesh()
  |-- 写出 transfer_function_data.h5
  |
  |-- AdjointFlux::compute_from_files()
  |-- 读取 transfer_function_data.h5 + fission_matrix.h5
  |-- 写出 adjoint_flux.h5
  |
  |-- FluxMesh::finalize()
  |-- 写出 flux_mesh.h5
  |
  |-- BetaEffective::compute_from_files()
  |-- 读取 flux_mesh.h5 + adjoint_flux.h5
  |-- 写出 beta_eff.h5
```

实际触发位置在 `src/simulation.cpp`：

- batch 初始化时创建 `FissionMatrix` 和 `FluxMesh`；
- `current_batch == n_inactive` 时计算并写出 `fission_matrix.h5`；
- `current_batch == n_batches && clutch_on` 时依次写出 transfer function、adjoint flux、flux mesh 和 beta_eff。

## 2. 正向通量输入：flux_mesh.h5

正向通量由 `FluxMesh` 在粒子移动时累积。

在 `src/particle.cpp` 的 `Particle::event_advance()` 中，只要 `settings::flux_mesh_on && simulation::flux_mesh`，每次粒子移动都会调用：

```cpp
simulation::flux_mesh->accumulate(
  position, wgt(), distance, energy_eV, mg_group);
```

也就是说，`flux_mesh.h5` 中的通量来自 track-length 累积。

`BetaEffective::read_flux_data()` 会从 `flux_mesh.h5` 读取：

- `grid_shape`
- `grid_pitch`
- `grid_lower_left`
- `n_groups`
- `energy_edges`
- `cell_indices`
- `flux_mean`
- `flux_group_mean`

其中：

- `flux_mean` 被放入 `flux_map[cell]`；
- `flux_group_mean` 被整理为 `flux_group_map_[cell][g]`；
- 当前 `validate_group_metadata()` 强制要求正向通量具备多群数据。

## 3. 伴随源输入：fission_matrix.h5

`FissionMatrix` 在 inactive batches 中运行。当前 `simulation.cpp` 里实际启用的是：

```cpp
if (simulation::fission_matrix &&
    simulation::current_batch <= settings::n_inactive) {
  simulation::fission_matrix->start_new_batch(simulation::current_batch);
}
```

inactive batches 结束时：

```cpp
simulation::fission_matrix->compute_adjoint_source(...);
simulation::fission_matrix->finalize("fission_matrix.h5");
```

`compute_adjoint_source()` 的关键过程是：

1. 用源计数构造经验源能谱 `chi_empirical(g|cell)`。
2. 将稀疏裂变矩阵按源状态计数归一化：

```text
M_norm[source_state -> fission_cell]
```

1. 初始化 `adjoint_source_grouped_`。
2. 做幂迭代：

```text
q(j) = (1 / keff) * sum_s M_norm[s -> j] * I*(s)
I*_new(j,g) = q(j) * chi_empirical(g|j)
```

1. 每轮归一化 `I*_new`，直到收敛或达到最大迭代数。
2. 将分群伴随源折叠成标量：

```text
adjoint_source_[cell] = sum_g adjoint_source_grouped_[cell,g]
```

`fission_matrix.h5` 后续给 `AdjointFlux` 使用，优先读取：

- `n_source_groups`
- `adjoint_source_grouped`

如果没有分群伴随源，则回退读取：

- `adjoint_source`

## 4. 传递函数输入：transfer_function_data.h5

传递函数由 `GreenFunctionMesh` 统计。

### 4.1 源状态记录

`GreenFunctionMesh::record_source_birth()` 将源粒子映射到源状态：

```text
source_state = source_cell * n_source_groups + g_source
```

该映射存入：

```text
particle_to_source_state_[source_particle_id] = source_state
```

后续同一源粒子的裂变贡献会通过 `source_particle_id` 追溯到这个源状态。

### 4.2 裂变事件贡献

在 `src/tallies/tally_scoring.cpp` 的 `SCORE_GREENFUNCTION` 分支中，只有发生裂变的中子才累积传递函数。

先计算总裂变贡献：

```cpp
double nu_t = p.wgt() / simulation::keff * weight *
              p.neutron_xs(p.event_nuclide()).nu_fission /
              p.neutron_xs(p.event_nuclide()).total;
```

然后按 family 拆分：

```text
family 0      = prompt
family 1..8   = delayed group 1..8
```

prompt 分量：

```cpp
nu_t * (nu_p / nu_total_yield)
```

delayed 第 `k` 群分量：

```cpp
nu_t * (nu_dk / nu_total_yield)
```

这些分量调用 `GreenFunctionMesh::accumulate()` 写入：

```text
T(source_state -> response_cell, family, group)
```

内部存储顺序是：

```text
flat_idx = family * n_groups + group
```

### 4.3 文件输出

最后一个 batch 结束时，`GreenFunctionMesh::finalize_greenfunction_mesh()` 写出 `transfer_function_data.h5`，核心字段包括：

- `origin`
- `shape`
- `pitch`
- `n_source_groups`
- `n_groups`
- `n_families`
- `family_order = prompt,delayed_1,...,delayed_8`
- `storage_order = family_major: [family][group]`
- `energy_edges`
- `source_state_indices`
- `source_state_cell`
- `source_state_group`
- `transfer_functions/source_state_N/indices`
- `transfer_functions/source_state_N/values`

每个 `values` 是按 response cell 存储的一段 family-major 向量。

## 5. 伴随通量文件：adjoint_flux.h5

`AdjointFlux::compute_from_files()` 读取：

```text
transfer_function_data.h5
fission_matrix.h5
```

并计算：

```text
Phi_dag(response_cell, family, group)
  = sum_source_state T(source_state -> response_cell, family, group)
                     * I*(source_state)
```

实现位置在 `AdjointFlux::compute_from_memory()`。

如果 `n_families > 1`，会同时累积：

- 总的 `adjoint_flux_sparse_[response_cell][group]`
- 分 family 的 `family_adjoint_flux_[family][response_cell][group]`

`adjoint_flux.h5` 写出：

- `cell_indices`
- `flux_mean`
- `flux_group_mean`
- `adjoint_flux_dense`
- `adjoint_flux_group_dense`
- `group_total_flux`
- `family_resolved/`

其中 `family_resolved/` 下的结构是：

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

注意：当前 `BetaEffective` 后续只读取 `family_resolved/*/flux_mean`，没有在 Method E 中使用 `family_resolved/*/flux_group_mean`。

## 6. BetaEffective 主流程

入口函数是：

```cpp
BetaEffective::compute_from_files(
  "flux_mesh.h5", "adjoint_flux.h5", "beta_eff.h5")
```

内部分四步：

1. 读取正向通量 `flux_mesh.h5`。
2. 读取伴随通量和 family-resolved importance `adjoint_flux.h5`。
3. 构建 cell 到材料核数据的映射。
4. 使用 Method E 计算分母、各群分子和 beta_eff，并写出 `beta_eff.h5`。

## 7. 输入检查

`BetaEffective::validate_group_metadata()` 当前强制要求：

- `flux_mesh.h5` 有多群正向通量；
- `adjoint_flux.h5` 有多群伴随通量；
- 两者 `n_groups` 一致；
- 两者能群边界 `energy_edges` 一致；
- `flux_group_map_` 和 `adjoint_group_map_` 非空。

否则直接 `fatal_error()`。

虽然 `compute_denominator_upstream_family()` 和 `compute_delayed_numerator_upstream_family()` 里保留了单群回退公式，但当前入口检查已经要求多群输入，因此正常路径是多群。

## 8. 材料核数据构建

`build_cell_material_map(flux)` 会根据正向通量中的非零 cell 建立材料映射，并为每个 cell 准备 `MaterialNuclearData`。

主要目的：

- 判断 cell 是否含裂变材料；
- 给每个 cell 提供：
  - `sigma_f_groups[g]`
  - `nu_total_groups[g]`
  - `nu_delayed_groups[g,k]`
  - 单群回退用的 `sigma_f`
  - 单群回退用的 `nu_total`
- 对异质 cell 使用多点采样和材料加权。

默认构造参数是：

```cpp
BetaEffective beta_calc;
```

即：

- `n_sample_points = 27`
- `weighting_mode = REACTION_RATE_WEIGHTED`
- `ref_energy = 2.0e6 eV`
- `ref_temperature = 293.6 K`

## 9. 当前主算法：Method E

当前 `beta_effective.cpp` 只把 Method E 作为主结果。

代码注释写明 Method E 是：

```text
上游族解析 upstream family-resolved adjoint
```

它要求 `adjoint_flux.h5` 中存在：

```text
family_resolved/
```

否则直接报错：

```text
方法E需要上游族解析伴随通量数据
```

### 9.1 Method E 读取的数据

`read_adjoint_flux_data()` 从 `family_resolved/` 读取：

prompt:

```text
family_resolved/prompt/cell_indices
family_resolved/prompt/flux_mean
```

delayed:

```text
family_resolved/delayed_1/cell_indices
family_resolved/delayed_1/flux_mean
...
family_resolved/delayed_8/cell_indices
family_resolved/delayed_8/flux_mean
```

整理到：

```cpp
cell_upstream_prompt_importance_[cell]
cell_upstream_delayed_importance_[k][cell]
```

### 9.2 分母计算

函数：

```cpp
compute_denominator_upstream_family(flux, volume)
```

对每个正向通量非零 cell：

1. 查 `cell_nuclear_data_`，跳过非裂变 cell。
2. 计算：

```text
I_total(cell) = I_prompt(cell) + sum_k I_delayed_k(cell)
```

1. 用多群正向通量和核数据计算：

```text
F_total(cell) = sum_g nu_total_g(cell) * sigma_f_g(cell) * phi_g(cell)
```

1. 累积：

```text
D += volume * I_total(cell) * F_total(cell)
```

最后用 `kahan_sum()` 返回分母。

### 9.3 分子计算

函数：

```cpp
compute_delayed_numerator_upstream_family(group, flux, volume)
```

对每个正向通量非零 cell：

1. 查 `cell_upstream_delayed_importance_[group][cell]`，没有则跳过。
2. 查 `cell_nuclear_data_`，跳过非裂变 cell。
3. 计算同一个：

```text
F_total(cell) = sum_g nu_total_g(cell) * sigma_f_g(cell) * phi_g(cell)
```

1. 累积：

```text
N_group += volume * I_delayed_group(cell) * F_total(cell)
```

注意：这里分子使用 `F_total`，不是 `F_delayed_group`。代码注释明确说明原因：`I_delayed_group` 已经在上游传递函数 family 轴中包含 `nu_d,k / nu_t` 产额分份，若再乘 `F_delayed_group` 会重复计数 delayed yield。

### 9.4 beta_i 和 beta_total

`compute_from_files()` 中循环 8 个 delayed group：

```cpp
num_upstream[i] = compute_delayed_numerator_upstream_family(i, flux, volume);
beta_upstream[i] = num_upstream[i] / denom_upstream;
beta_upstream_total += beta_upstream[i];
```

然后把 Method E 结果作为主结果：

```cpp
denominator_ = denom_upstream;
numerators_[i] = num_upstream[i];
beta_i_[i] = beta_upstream[i];
beta_total_ = beta_upstream_total;
```

如果 `denominator_ <= 0.0`，直接报错。

## 10. 输出文件：beta_eff.h5

`BetaEffective::write_to_file()` 输出 `beta_eff.h5`。

文件属性：

- `filetype = beta_effective`
- `version = 1.0`
- `description = Effective delayed neutron fraction (β_eff)`

主结果：

```text
beta_i              # 8 组 beta_i
beta_total          # 总 beta_eff
uncertainty         # 当前写 8 个 0
```

metadata:

```text
metadata/flux_file
metadata/adjoint_flux_file
metadata/computation_time
metadata/nuclear_data_source
metadata/energy_groups
metadata/energy_edges
metadata/beta_eff_mode = MATERIAL_DEPENDENT
metadata/ref_energy_ev
metadata/ref_temperature_k
metadata/n_delayed_groups = 8
metadata/grid_shape
metadata/grid_pitch
```

diagnostics:

```text
diagnostics/numerator
diagnostics/denominator
diagnostics/normalization_factor
diagnostics/material_ids
diagnostics/material_nu_total
diagnostics/material_sigma_f
diagnostics/material_cell_counts
diagnostics/n_unique_materials
diagnostics/total_fissionable_cells
```

method comparison:

```text
method_comparison/method_e/beta_i
method_comparison/method_e/beta_total
method_comparison/method_e/numerator
method_comparison/method_e/denominator
```

当前只有 `method_e` 被写入。

## 11. 当前代码流程中的几个重要事实

1. `beta_eff` 当前依赖 `clutch_on` 链路生成的三个文件：

```text
flux_mesh.h5
adjoint_flux.h5
fission_matrix.h5
```

其中 `fission_matrix.h5` 不被 `BetaEffective` 直接读取，但它是生成 `adjoint_flux.h5` 的必要输入。

1. `BetaEffective` 入口要求 `flux_mesh.h5` 和 `adjoint_flux.h5` 都有多群数据。

2. `BetaEffective` 的 Method E 使用的是 `family_resolved/*/flux_mean` 标量 family importance。

3. `adjoint_flux.h5` 虽然写出了 `family_resolved/*/flux_group_mean`，但当前 Method E 没有使用这些 family 内部分群数据。

4. `cell_upstream_delayed_importance_[k]` 已经包含上游按 `nu_d,k / nu_total` 拆分后的 delayed family 贡献，所以分子乘的是 `F_total`。

5. `uncertainty` 当前固定写零，代码中还没有对 `beta_eff` 做 batch-wise ratio estimator 不确定度传播。

6. `beta_effective.cpp` 中还打印了 MCNP 参考值对比，但这些参考值只用于控制台对比，不参与 `beta_eff` 计算。

## 12. 用一句话概括当前实现

当前代码中的 `beta_eff` 计算是：

```text
先用 fission_matrix.h5 的 I*(source_state)
与 transfer_function_data.h5 的 family-resolved T(source_state -> cell)
卷积得到 adjoint_flux.h5 中的 prompt/delayed family importance；
再用 flux_mesh.h5 的多群正向通量和每个 cell 的材料核数据构造 F_total；
最后按
  D_k? no: D = sum_cell volume * (I_prompt + sum I_delayed_k) * F_total
  N_k     = sum_cell volume * I_delayed_k * F_total
  beta_k  = N_k / D
求得 8 组 beta_i 和 beta_total。
```
