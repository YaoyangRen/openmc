# 广义 CLUTCH 敏感性系数 XML 输入使用手册

本文说明如何通过 XML 输入启用当前实现的广义 CLUTCH k-effective 敏感性系数计算。

## 功能范围

当前首版实现用于计算 `k-effective` 对材料扰动参数的 CLUTCH 敏感性系数。

支持：

- 连续能量 `CE` 本征值计算。
- history-based 输运。
- `openmc.TallyDerivative` 已支持的三类扰动：
  - `density`
  - `nuclide_density`
  - `temperature`
- 主方法 `fclutch_fm`：使用 inactive batches 形成的裂变矩阵空间伴随源 `I*(cell)`。
- 第二方法 `cclutch_history`：空间 transfer-function C-CLUTCH，对 active batch 中的源 cell 到响应 cell 裂变生产传递函数进行折叠。

暂不支持：

- MG 模式。
- event-based 模式。
- fixed-source 模式。
- MT/能群级核数据敏感性。
- 出生能量维度的 C-CLUTCH。

## 必需输入文件

启用该功能至少需要：

- `settings.xml` 中新增 `<clutch_sensitivity>`。
- `tallies.xml` 中定义一个或多个 `<derivative>`。
- `settings.xml` 中有 inactive batches，且 `inactive > 0`。

`<clutch_sensitivity>` 出现后，程序会自动启用裂变矩阵 CLUTCH 源重要性计算。无需添加 `clutch-test` 或 `greenfunction` tally。

## settings.xml 示例

最小示例：

```xml
<settings>
  <run_mode>eigenvalue</run_mode>
  <particles>10000</particles>
  <batches>120</batches>
  <inactive>40</inactive>

  <kinetics_mesh>
    <pitch>1.0</pitch>
    <auto_bounds>false</auto_bounds>
    <lower_left>-10.0 -10.0 -10.0</lower_left>
    <upper_right>10.0 10.0 10.0</upper_right>
  </kinetics_mesh>

  <adjoint_source>
    <initial_guess>uniform</initial_guess>
    <max_iterations>100</max_iterations>
    <tolerance>1.0e-8</tolerance>
    <score_start_batch>3</score_start_batch>
  </adjoint_source>

  <clutch_sensitivity>
    <method>hybrid</method>
    <output>clutch_sensitivity.h5</output>
    <derivative_ids>1 2 3</derivative_ids>
  </clutch_sensitivity>
</settings>
```

## clutch_sensitivity 节点

`<clutch_sensitivity>` 是启用节点。只要该节点存在，就会启用广义 CLUTCH sensitivity。

字段：

| 字段 | 是否必需 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `method` | 否 | `hybrid` | 输出方法，可选 `hybrid`、`fclutch_fm`、`cclutch_history`。 |
| `output` | 否 | `clutch_sensitivity.h5` | sensitivity HDF5 输出文件名。 |
| `derivative_ids` | 否 | 使用所有 `<derivative>` | 指定参与计算的 derivative ID 列表。 |

`method` 说明：

- `hybrid`：同时写出 `/method/fclutch_fm` 和 `/method/cclutch_history`。
- `fclutch_fm`：只写主方法结果。
- `cclutch_history`：只写空间 transfer-function C-CLUTCH 结果。

建议首选 `hybrid`，因为它能同时看到 F-CLUTCH 和 C-CLUTCH 两种方法的物理对照。

## kinetics_mesh 节点

`<kinetics_mesh>` 定义裂变矩阵和空间伴随源 `I*(cell)` 的空间网格。

字段：

| 字段 | 是否必需 | 说明 |
| --- | --- | --- |
| `pitch` | 否 | 网格尺寸，单位 cm，默认 `1.0`。 |
| `auto_bounds` | 否 | 是否自动从几何边界取范围，默认 `true`。 |
| `lower_left` | 条件必需 | 当 `auto_bounds=false` 时必须给出。 |
| `upper_right` | 条件必需 | 当 `auto_bounds=false` 时必须给出。 |

建议：

- 对调试和有限差分比较，先用较粗网格，保证 inactive 阶段每个重要裂变区域都有足够采样。
- 若几何有无限边界，必须手动给出 `lower_left` 和 `upper_right`。

示例：

```xml
<kinetics_mesh>
  <pitch>0.5</pitch>
  <auto_bounds>false</auto_bounds>
  <lower_left>-5.0 -5.0 -1.0</lower_left>
  <upper_right>5.0 5.0 1.0</upper_right>
</kinetics_mesh>
```

## adjoint_source 节点

`<adjoint_source>` 控制 inactive batches 裂变矩阵伴随源迭代。

字段：

| 字段 | 是否必需 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `initial_guess` | 否 | `uniform` | 初始猜测，可选 `uniform`、`forward`。 |
| `max_iterations` | 否 | `1` | 伴随源迭代最大步数。实际 sensitivity 建议设置为 `50-200`。 |
| `tolerance` | 否 | `1.0e-6` | 伴随源迭代收敛阈值。 |
| `score_start_batch` | 否 | `3` | 从第几个 inactive batch 开始累计裂变矩阵。 |

示例：

```xml
<adjoint_source>
  <initial_guess>forward</initial_guess>
  <max_iterations>100</max_iterations>
  <tolerance>1.0e-8</tolerance>
  <score_start_batch>5</score_start_batch>
</adjoint_source>
```

`score_start_batch` 用于跳过 inactive 初期尚未收敛的源分布。若设置值大于 `inactive`，程序会自动降到 `inactive` 并给出 warning。

## tallies.xml 中定义 derivative

CLUTCH sensitivity 使用 `openmc.TallyDerivative` 的 XML 定义，但不要求普通 tally 引用这些 derivative。

示例：

```xml
<tallies>
  <derivative id="1" variable="density" material="1" />
  <derivative id="2" variable="nuclide_density" material="1" nuclide="U235" />
  <derivative id="3" variable="temperature" material="1" />
</tallies>
```

字段：

| 属性 | 是否必需 | 说明 |
| --- | --- | --- |
| `id` | 是 | derivative ID。用于 `<derivative_ids>` 选择。 |
| `variable` | 是 | `density`、`nuclide_density` 或 `temperature`。 |
| `material` | 是 | 被扰动材料 ID。 |
| `nuclide` | 条件必需 | 仅 `variable="nuclide_density"` 时需要。 |

如果 `<clutch_sensitivity>` 不写 `<derivative_ids>`，程序默认使用 `tallies.xml` 中所有 `<derivative>`。

## 完整 XML 示例

`settings.xml`：

```xml
<?xml version="1.0"?>
<settings>
  <run_mode>eigenvalue</run_mode>
  <particles>20000</particles>
  <batches>150</batches>
  <inactive>50</inactive>

  <energy_mode>continuous-energy</energy_mode>
  <event_based>false</event_based>

  <kinetics_mesh>
    <pitch>1.0</pitch>
    <auto_bounds>false</auto_bounds>
    <lower_left>-10.0 -10.0 -10.0</lower_left>
    <upper_right>10.0 10.0 10.0</upper_right>
  </kinetics_mesh>

  <adjoint_source>
    <initial_guess>uniform</initial_guess>
    <max_iterations>100</max_iterations>
    <tolerance>1.0e-8</tolerance>
    <score_start_batch>5</score_start_batch>
  </adjoint_source>

  <clutch_sensitivity>
    <method>hybrid</method>
    <output>clutch_sensitivity.h5</output>
    <derivative_ids>1 2</derivative_ids>
  </clutch_sensitivity>
</settings>
```

`tallies.xml`：

```xml
<?xml version="1.0"?>
<tallies>
  <derivative id="1" variable="density" material="1" />
  <derivative id="2" variable="nuclide_density" material="1" nuclide="U235" />
</tallies>
```

## Python 生成 XML 示例

推荐通过 Python API 生成 XML，避免手写 ID 或节点格式错误。

```python
import openmc

model = openmc.Model()

density_deriv = openmc.TallyDerivative(
    derivative_id=1,
    variable='density',
    material=1
)

u235_deriv = openmc.TallyDerivative(
    derivative_id=2,
    variable='nuclide_density',
    material=1,
    nuclide='U235'
)

model.settings.particles = 20000
model.settings.batches = 150
model.settings.inactive = 50
model.settings.energy_mode = 'continuous-energy'
model.settings.event_based = False

model.enable_clutch_sensitivity(
    derivatives=[density_deriv, u235_deriv],
    method='hybrid',
    output='clutch_sensitivity.h5',
    kinetics_mesh={
        'pitch': 1.0,
        'auto_bounds': False,
        'lower_left': (-10.0, -10.0, -10.0),
        'upper_right': (10.0, 10.0, 10.0),
    },
    adjoint_source={
        'initial_guess': 'uniform',
        'max_iterations': 100,
        'tolerance': 1.0e-8,
        'score_start_batch': 5,
    }
)

model.export_to_xml()
```

注意：上例只展示 sensitivity 相关设置。实际算例仍需正常定义 `materials.xml`、`geometry.xml`、source 等输入。

## 输出文件结构

默认输出文件为 `clutch_sensitivity.h5`。

主要内容：

```text
/parameters
  ids
  variable
  material_id
  nuclide
  parameter_value

/method/fclutch_fm
  dlogk_dparameter
  sensitivity
  numerator
  denominator
  uncertainty
  batch_ids
  batch_numerator
  batch_denominator

/method/cclutch_history
  dlogk_dparameter
  sensitivity
  numerator
  denominator
  uncertainty
  batch_ids
  batch_numerator
  batch_denominator

/diagnostics
  grid_shape
  grid_lower_left
  grid_upper_right
  grid_pitch
```

物理含义：

- `dlogk_dparameter`：`d ln(k) / d parameter`。
- `sensitivity`：无量纲敏感性，等于 `parameter_value * dlogk_dparameter`。
- `uncertainty`：对 `sensitivity` 的 batch 统计不确定度。
- `fclutch_fm` 的分母来自 active fission source site 的 `w_site * I*(cell_site)`。
- `cclutch_history` 的分母来自 active transfer function：`sum_source I*(source_cell) * T_total(source_cell)`。

## 分析脚本输出解读

可以使用仓库根目录下的脚本分析输出文件：

```bash
python py_script/analyze_clutch_sensitivity.py clutch_sensitivity.h5
```

脚本输出中的 `Method: fclutch_fm` 是主方法结果，`Method: cclutch_history (transfer-function C-CLUTCH)` 是第二方法结果。两者使用同一个 inactive 裂变矩阵伴随源 `I*(cell)`，但 active 阶段的计分对象不同：F-CLUTCH 直接在裂变源 site 上加权，C-CLUTCH 先统计源 cell 到响应 cell 的裂变生产 transfer function，再用 `I*(source_cell)` 折叠。

### 主结果表

示例：

```text
   id         variable    mat      nuclide      parameter       dlogk/dp    sensitivity          unc      rel_unc
    1          density      1            -   1.875000e+01   4.409819e-02   8.268411e-01   1.7376e-03    2.102e-03
    2  nuclide_density      1         U235   4.501828e-02   1.779730e+01   8.012037e-01   1.7302e-03    2.160e-03
    3      temperature      1            -   2.936000e+02   0.000000e+00   0.000000e+00   0.0000e+00          nan
```

字段解释：

- `id`：`tallies.xml` 中 `<derivative>` 的 ID。
- `variable`：扰动类型，包括 `density`、`nuclide_density`、`temperature`。
- `mat`：被扰动材料 ID。
- `nuclide`：被扰动核素。`density` 和 `temperature` 对整个材料扰动，因此显示为 `-`。
- `parameter`：扰动参数当前值。例如材料密度、核素原子密度、材料温度。
- `dlogk/dp`：`d ln(k) / d parameter`，带有参数单位的导数。
- `sensitivity`：无量纲敏感性，等于 `parameter * dlogk/dp`。
- `unc`：`sensitivity` 的统计不确定度。
- `rel_unc`：相对不确定度，等于 `abs(unc / sensitivity)`。

若输出为：

```text
density sensitivity = 8.268411e-01
U235 nuclide_density sensitivity = 8.012037e-01
```

可近似理解为：

- 材料 1 总密度增加 1%，`k_eff` 增加约 `0.8268%`。
- 材料 1 中 `U235` 核素密度增加 1%，`k_eff` 增加约 `0.8012%`。

### denominator

示例：

```text
denominator mean: 3.78916978e+00
```

这是 CLUTCH ratio estimator 的归一化分母。对 `fclutch_fm`：

```text
D = sum_fission_site w_site * I*(cell_site)
```

对 `cclutch_history`：

```text
D = sum_source_cell I*(source_cell) * T_total(source_cell)
```

其中 `T_total(source_cell)` 是该 active batch 中从该源 cell 出发的单位源粒子平均总裂变生产响应。它不是敏感性本身，而是归一化响应量。分母越稳定，说明 active source/transfer function 与 inactive 生成的空间重要性匹配越稳定。

### Batch diagnostics

示例：

```text
batches: 200
zero denominator batches: 0
denominator min/mean/max: 3.656722e+00 / 3.789170e+00 / 3.932452e+00
denominator coeff. variation: 1.515877e-02
```

含义：

- `batches`：参与统计的 active batch 数。
- `zero denominator batches`：分母为零的 batch 数。正常应为 `0`。
- `denominator min/mean/max`：每个 batch 的分母范围。
- `denominator coeff. variation`：分母批间变异系数，越小表示归一化响应越稳定。

如果 denominator 变异系数很大，通常说明 inactive 裂变矩阵伴随源采样不足、网格过细，或 active fission source 与伴随源非零区域匹配不好。

### Batch ratio

示例：

```text
batch ratio dlogk/dp mean/std by parameter:
  param[01]: mean=4.409295e-02, std=1.314935e-03
  param[02]: mean=1.779491e+01, std=5.452429e-01
```

这里显示每个 batch 的 `N_b / D_b` 的均值和标准差，用于观察批间波动。最终主表中的 `dlogk/dp` 使用 ratio-of-means：

```text
mean(N_b) / mean(D_b)
```

因此它通常与 batch ratio mean 很接近，但不要求完全相等。

### temperature sensitivity 为零

若看到：

```text
temperature dlogk/dp = 0
temperature sensitivity = 0
```

通常表示当前 temperature derivative 没有可用贡献。常见原因：

- 材料相关核素没有 Windowed Multipole 数据。
- 粒子能量不在 multipole 支持能区内。
- 未启用或未加载 temperature derivative 所需的 multipole 数据。

当前实现中，`temperature` sensitivity 依赖 multipole 数据计算 `d sigma / dT`。如果核素没有 multipole 数据，程序会将对应 temperature 贡献按 0 处理。

### Method comparison

示例：

```text
Method comparison: C-CLUTCH - F-CLUTCH
id  variable          delta_sens      rel_delta
1   density           0.000000e+00    0.000000e+00
2   nuclide_density   0.000000e+00    0.000000e+00
```

现在 `cclutch_history` 是真正的空间 transfer-function C-CLUTCH，因此 `delta_sens` 表示 C-CLUTCH 与 F-CLUTCH 两种估计路线的差值。若差值明显大于统计不确定度，优先检查 kinetics mesh 是否过细、inactive 伴随源覆盖率是否不足、以及 C-CLUTCH denominator 的 batch 变异系数是否偏大。

## 常见错误

### 缺少 derivative

如果启用了 `<clutch_sensitivity>`，但 `tallies.xml` 中没有 `<derivative>`，会报错：

```text
CLUTCH sensitivity requires at least one <derivative> in tallies.xml.
```

解决：添加至少一个 `<derivative>`，或用 `Model.enable_clutch_sensitivity()` 自动生成。

### derivative_ids 指向不存在 ID

如果 `<derivative_ids>` 中的 ID 没有在 `tallies.xml` 定义，会报错。

解决：检查 `<derivative id="...">` 与 `<derivative_ids>` 是否一致。

### event-based / MG / fixed-source

当前实现会直接拒绝这些模式。

解决：使用 CE eigenvalue history-based：

```xml
<run_mode>eigenvalue</run_mode>
<energy_mode>continuous-energy</energy_mode>
<event_based>false</event_based>
```

### inactive 为 0

该方法需要 inactive batches 构建裂变矩阵和空间伴随源。

解决：设置：

```xml
<inactive>至少为 1</inactive>
```

实际计算建议 inactive 足够大，例如 `30-100`，视算例规模而定。

### sensitivity denominator 为 0

通常表示 active fission sites 没有落在非零空间伴随源区域，或网格过细导致 inactive 裂变矩阵欠采样。

建议：

- 增加 inactive/active histories。
- 粗化 `kinetics_mesh/pitch`。
- 检查 `kinetics_mesh` 边界是否覆盖所有裂变区域。
- 查看 `fission_matrix.h5` 和 `clutch_sensitivity.h5:/diagnostics`。

## 推荐工作流程

1. 先用较粗 `kinetics_mesh` 跑小统计，确认能生成非零 `fission_matrix.h5` 和 `clutch_sensitivity.h5`。
2. 查看 `/diagnostics/total_scored_sites` 和 `/diagnostics/total_dropped_sites`。
3. 若 dropped sites 较多，检查网格边界。
4. 若不确定度大，增加 active batches 或 histories。
5. 用有限差分扰动材料密度或核素密度，对比 `sensitivity` 是否在统计不确定度内一致。
