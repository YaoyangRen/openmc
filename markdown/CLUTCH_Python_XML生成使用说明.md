# CLUTCH Python XML 生成使用说明

## 目的

本文档说明如何使用 OpenMC 的 Python API 生成 CLUTCH/beta-effective 计算所需的 `settings.xml` 和 `tallies.xml`。目标是避免手写 XML，只需要在 Python 脚本中设置 `clutch` 开启以及相关参数，即可生成与当前手写输入等价的 XML 文件。

当前实现采用以下约定：

- `settings.xml` 中写入动力学计算参数，例如 `kinetics_energy_edges`、`kinetics_mesh` 和 `adjoint_source`。
- `tallies.xml` 中写入 CLUTCH 所需 score，例如 `clutch-test` 和 `greenfunction`。
- C++ 侧没有单独读取 `<clutch_on>true</clutch_on>` 节点；当 tally 中出现 `clutch-test` 或 `greenfunction` 时，输入读取阶段会自动设置 `settings::clutch_on = true`，并同时启用 `flux_mesh_on`。

因此，Python API 中的一键接口 `model.enable_clutch(...)` 实际完成两件事：

1. 设置 `Settings` 中的动力学参数。
2. 向 `Tallies` 中添加开启 CLUTCH 所需的 tally score。

## 推荐用法：Model.enable_clutch

推荐使用 `openmc.Model.enable_clutch(...)`。该接口会同时配置 `settings.xml` 和 `tallies.xml`。

```python
import openmc
import openmc.stats

model = openmc.Model()

model.settings.run_mode = 'eigenvalue'
model.settings.particles = 1000
model.settings.batches = 255
model.settings.inactive = 55
model.settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0.0, 0.0, 0.0))
)

energy_edges = [
    1.000010E-05, 1.000000E-01, 5.400000E-01, 4.000000E+00,
    8.315287E+00, 1.370959E+01, 2.260329E+01, 4.016900E+01,
    6.790405E+01, 9.166088E+01, 1.486254E+02, 3.043248E+02,
    4.539993E+02, 7.485183E+02, 1.234098E+03, 2.034684E+03,
    3.354626E+03, 5.530844E+03, 9.118820E+03, 1.503439E+04,
    2.478752E+04, 4.086771E+04, 6.737947E+04, 1.110900E+05,
    1.831564E+05, 3.019738E+05, 4.978707E+05, 8.208500E+05,
    1.353353E+06, 2.231302E+06, 3.678794E+06, 6.065307E+06,
    1.000000E+07, 1.964033E+07,
]

model.enable_clutch(
    kinetics_energy_edges=energy_edges,
    adjoint_source={
        'initial_guess': 'uniform',
        'max_iterations': 100,
        'tolerance': 1.0e-8,
    }
)

model.settings.export_to_xml('build/bin/beta_test_inp')
model.tallies.export_to_xml('build/bin/beta_test_inp')
```

该脚本会生成：

```text
build/bin/beta_test_inp/settings.xml
build/bin/beta_test_inp/tallies.xml
```

## 显式指定动力学网格

如果需要指定 `<kinetics_mesh>`，可在 `enable_clutch()` 中传入 `kinetics_mesh` 字典：

```python
model.enable_clutch(
    kinetics_energy_edges=energy_edges,
    kinetics_mesh={
        'pitch': 1.0,
        'auto_bounds': True,
    },
    adjoint_source={
        'initial_guess': 'uniform',
        'max_iterations': 100,
        'tolerance': 1.0e-8,
    }
)
```

若几何边界为无限，或者希望完全手动指定网格范围：

```python
model.enable_clutch(
    kinetics_energy_edges=energy_edges,
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
    }
)
```

对应 XML：

```xml
<kinetics_mesh>
  <pitch>1.0</pitch>
  <auto_bounds>false</auto_bounds>
  <lower_left>-10.0 -10.0 -10.0</lower_left>
  <upper_right>10.0 10.0 10.0</upper_right>
</kinetics_mesh>
```

## 只使用 CLUTCH score，不写 IFP score

默认情况下，`enable_clutch()` 会生成当前测试输入中的完整 score 列表：

```text
ifp-time-numerator ifp-beta-numerator ifp-denominator clutch-test greenfunction
```

如果只想生成 CLUTCH 相关 score：

```python
model.enable_clutch(
    kinetics_energy_edges=energy_edges,
    include_ifp=False,
)
```

此时 `tallies.xml` 中只包含：

```text
clutch-test greenfunction
```

## 生成的 tallies.xml

默认调用：

```python
model.enable_clutch(kinetics_energy_edges=energy_edges)
```

会生成类似下面的 `tallies.xml`：

```xml
<?xml version='1.0' encoding='utf-8'?>
<tallies>
  <tally id="1" name="scores">
    <scores>ifp-time-numerator ifp-beta-numerator ifp-denominator clutch-test greenfunction</scores>
  </tally>
</tallies>
```

其中 `greenfunction` 用于传递函数统计，`clutch-test` 用于触发 CLUTCH 相关逻辑。C++ 输入读取器看到这两个 score 后，会自动启用：

```text
settings::clutch_on = true
settings::flux_mesh_on = true
```

## 生成的 settings.xml

典型 `settings.xml` 会包含：

```xml
<settings>
  <run_mode>eigenvalue</run_mode>
  <particles>1000</particles>
  <batches>255</batches>
  <inactive>55</inactive>
  <source strength="1.0">
    <space type="point">
      <parameters>0.0 0.0 0.0</parameters>
    </space>
  </source>
  <kinetics_energy_edges>...</kinetics_energy_edges>
  <adjoint_source>
    <initial_guess>uniform</initial_guess>
    <max_iterations>100</max_iterations>
    <tolerance>1e-08</tolerance>
  </adjoint_source>
</settings>
```

如传入 `kinetics_mesh`，还会额外生成：

```xml
<kinetics_mesh>
  ...
</kinetics_mesh>
```

## 直接使用 Settings 和 Tallies

如果不使用 `Model.enable_clutch()`，也可以手动调用底层 API：

```python
import openmc
import openmc.stats

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.particles = 1000
settings.batches = 255
settings.inactive = 55
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0.0, 0.0, 0.0))
)
settings.kinetics_energy_edges = energy_edges
settings.adjoint_source = {
    'initial_guess': 'uniform',
    'max_iterations': 100,
    'tolerance': 1.0e-8,
}
settings.kinetics_mesh = {
    'pitch': 1.0,
    'auto_bounds': True,
}

tallies = openmc.Tallies()
tallies.add_clutch_tally()

settings.export_to_xml('build/bin/beta_test_inp')
tallies.export_to_xml('build/bin/beta_test_inp')
```

## 参数说明

### kinetics_energy_edges

类型：

```python
Iterable[float]
```

要求：

- 至少两个边界。
- 所有值必须大于 0。
- 必须严格递增。

用于：

- `FissionMatrix` 的源状态能群划分。
- `GreenFunctionMesh` 的源能群和响应能群划分。
- `FluxMesh` 的分群通量统计。
- `BetaEffective` 的多群 `beta_eff` 计算。

### kinetics_mesh

类型：

```python
dict
```

允许 key：

```text
pitch
auto_bounds
lower_left
upper_right
```

建议：

- 有限几何优先使用 `auto_bounds=True`。
- 无限几何或需要严格控制范围时使用 `auto_bounds=False` 并提供 `lower_left/upper_right`。

### adjoint_source

类型：

```python
dict
```

允许 key：

```text
initial_guess
max_iterations
tolerance
```

其中：

- `initial_guess` 可为 `'uniform'` 或 `'forward'`。
- `max_iterations` 必须为正整数。
- `tolerance` 必须为正实数。

如果 `enable_clutch()` 中不显式传入 `adjoint_source`，默认值为：

```python
{
    'initial_guess': 'uniform',
    'max_iterations': 100,
    'tolerance': 1.0e-8,
}
```

## 与当前手写 XML 的对应关系

当前手写 `settings.xml` 中的：

```xml
<kinetics_energy_edges>...</kinetics_energy_edges>
<adjoint_source>
  <initial_guess>uniform</initial_guess>
  <max_iterations>100</max_iterations>
  <tolerance>1.0e-8</tolerance>
</adjoint_source>
```

对应 Python：

```python
model.enable_clutch(
    kinetics_energy_edges=energy_edges,
    adjoint_source={
        'initial_guess': 'uniform',
        'max_iterations': 100,
        'tolerance': 1.0e-8,
    }
)
```

当前手写 `tallies.xml` 中的：

```xml
<tally id="1" name="scores">
  <scores>ifp-time-numerator ifp-beta-numerator ifp-denominator clutch-test greenfunction</scores>
</tally>
```

对应 Python：

```python
model.tallies.add_clutch_tally()
```

或者由 `model.enable_clutch(...)` 自动完成。

## 常见注意事项

1. 不需要手写 `<clutch_on>true</clutch_on>`。当前 C++ 逻辑通过 tally score 自动开启 CLUTCH。
2. 如果只导出 `settings.xml` 而没有导出 `tallies.xml`，CLUTCH 不会被自动开启。
3. `kinetics_energy_edges` 必须和后续 beta 计算所需的能群结构一致。
4. 若使用 `auto_bounds=True` 但几何 bounding box 无限，需要在 `kinetics_mesh` 中提供手动 bounds，或设置 `auto_bounds=False`。
5. 生成 XML 后可直接运行 OpenMC；程序会在初始化阶段输出动力学网格的来源和划分结果。
