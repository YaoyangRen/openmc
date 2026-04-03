# 通量网格问题修复说明

## 问题诊断

**症状**: flux_mesh.h5 文件被创建,但 cell_indices 和 flux_mean 数组为空(长度为0)

**根本原因**:
通量累积代码原本放在 `score_tracklength_tally()` 函数中,但该函数只有在存在 tracklength 类型的 tally 时才会被调用:

```cpp
// src/particle.cpp
if (!model::active_tracklength_tallies.empty()) {
    score_tracklength_tally(*this, distance);
}
```

如果用户没有定义任何 tracklength tally,`score_tracklength_tally()` 永远不会执行,flux_mesh 的 `accumulate()` 也就永远不会被调用。

## 解决方案

将 flux_mesh 累积代码从 `tally_scoring.cpp::score_tracklength_tally()` 移到 `particle.cpp` 中粒子移动的主循环中,确保每次粒子移动时都会累积通量,无论是否有 tally。

### 修改的文件

1. **src/particle.cpp**
   - 添加头文件: `#include "openmc/flux_mesh.h"`
   - 在粒子移动后立即累积通量(在 tally 评分之前)

2. **src/tallies/tally_scoring.cpp**
   - 移除 flux_mesh 累积代码(避免重复)
   - 移除 `#include "openmc/flux_mesh.h"` (不再需要)

### 关键代码变更

**src/particle.cpp** (在 event_advance 函数中):

```cpp
// Advance particle in space and time
this->move_distance(distance);
double dt = distance / speed;
this->time() += dt;
this->lifetime() += dt;

// Accumulate flux mesh (if enabled) - 在所有粒子移动时累积
if (settings::flux_mesh_on && simulation::flux_mesh) {
  std::array<double, 3> position = {r().x, r().y, r().z};
  simulation::flux_mesh->accumulate(position, wgt(), distance);
}

// Score track-length tallies
if (!model::active_tracklength_tallies.empty()) {
  score_tracklength_tally(*this, distance);
}
```

## 测试步骤

1. 重新编译:

```bash
cd build
mingw32-make -j24
```

2. 运行模拟(确保 clutch_on = true):

```bash
./bin/openmc
```

3. 检查 flux_mesh.h5:

```python
import h5py
f = h5py.File('flux_mesh.h5', 'r')
print('非零单元数:', len(f['cell_indices'][:]))
print('通量均值样本:', f['flux_mean'][:5])
f.close()
```

4. 可视化:

```bash
python visualize_flux_mesh.py
```

## 预期结果

- flux_mesh.h5 应包含非零的 cell_indices 和 flux_mean 数据
- 可视化脚本应成功生成通量分布图
- 控制台应显示通量统计信息(最大值、最小值、平均值等)

## 技术细节

**为什么这样修改是正确的**:

1. **独立性**: flux_mesh 统计不应依赖于用户是否定义了 tally
2. **完整性**: 每次粒子移动都应该被统计(径迹长度估计器)
3. **性能**: 累积操作非常轻量,对性能影响可忽略
4. **一致性**: 与 fission_matrix 和 greenfunction_mesh 的处理方式一致

**位置选择**:

- ✅ 放在 `particle.cpp::event_advance()` 中粒子移动后
- ❌ ~~放在 `tally_scoring.cpp::score_tracklength_tally()` 中~~ (依赖 tally 存在)
- ❌ ~~放在其他 tally 函数中~~ (依赖特定 tally 类型)

---

**修复日期**: 2024-11-24  
**影响**: 所有使用 flux_mesh 功能的模拟
