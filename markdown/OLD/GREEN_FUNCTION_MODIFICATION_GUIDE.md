# 传递函数（Transfer Function）使用指南

**注意**: 本文档已更新以反映当前的实现。原先的"双格林函数网格"已简化为单一的"传递函数"。

---

## 📊 传递函数的物理意义

### 传递函数 (Transfer Function)

- **物理意义**: 从初始相空间点 P₀ 出发的源中子在空间位置 r 处产生的平均裂变中子数的期望值
- **数学定义**: `T(P₀ → r) = ∫∫ ν̄Σf(r,E')Φ(r,Ω',E'|S₀) dΩ'dE'`
- **计算公式**: `contribution = nu_t = (w/k_eff) × w_ufs × (ν̄Σf/Σt)`
- **累积时机**: 发生裂变事件时（在 `tally_scoring.cpp` 的 `SCORE_GREENFUNCTION` case 中）
- **应用**:
  - 伴随通量计算
  - 重要性函数求解
  - 扰动分析
  - 裂变链反应分析

### 实现特点

- ✅ **统计期望裂变中子数**: 使用 `nu_t`（期望值）而非采样后的整数
- ✅ **与裂变矩阵一致**: 使用相同的物理量和计算位置
- ✅ **Tally 控制开关**: 通过定义包含 `greenfunction` score 的 tally 来启用统计
- ✅ **代码维护性好**: 统计逻辑集中在 tally 系统中

---

## 🔧 使用方法

### 步骤 1: 在输入文件中定义 Tally

传递函数通过 tally 系统控制。在 `tallies.xml` 中添加：

```xml
<tallies>
  <tally id="1">
    <!-- 可选：添加过滤器 -->
    <filter type="cell" bins="1"/>
    
    <!-- 关键：包含 greenfunction score -->
    <scores>greenfunction</scores>
  </tally>
</tallies>
```

**说明**:

- 如果定义了包含 `greenfunction` score 的 tally，传递函数统计将自动启用
- 如果没有定义此 tally，传递函数将不会被计算，节省计算资源
- tally 本身返回 `score = 0.0`，只作为开关使用

### 步骤 2: 初始化传递函数网格

在 `src/simulation.cpp` 的 `initialize_batch()` 中：

```cpp
// Initialize transfer function mesh (传递函数)
if (!simulation::transfer_function_mesh) {
  simulation::transfer_function_mesh = std::make_unique<GreenFunctionMesh>(
    1.0,                    // 1cm 网格分辨率
    settings::n_batches,    // 最大批次数
    true                    // 自动获取边界
  );
}

```cpp
// 原有代码
std::unique_ptr<GreenFunctionMesh> green_function_mesh;

// 新增代码
std::unique_ptr<GreenFunctionMesh> fission_green_function_mesh;
```

#### 2.2 修改初始化代码

**位置**: `initialize_batch()` 函数中，约第 396-408 行

**原代码**:

```cpp
// Intialize greenfunction mesh
if (!simulation::green_function_mesh) {
  simulation::green_function_mesh = std::make_unique<GreenFunctionMesh>(
    1.0, settings::n_batches); // 1cm分辨率
}

// Green function 计数
if (settings::clutch_on && simulation::green_function_mesh) {
  simulation::green_function_mesh->start_new_batch(simulation::current_batch);
}
```

**修改为**:

```cpp
// Initialize flux green function mesh (通量格林函数)
if (!simulation::green_function_mesh) {
  simulation::green_function_mesh = std::make_unique<GreenFunctionMesh>(
    1.0,                    // 1cm 网格分辨率
    settings::n_batches,    // 最大批次数
    true                    // 自动获取边界
  );
}

// Initialize fission source green function mesh (裂变源格林函数)
if (!simulation::fission_green_function_mesh) {
  simulation::fission_green_function_mesh = std::make_unique<GreenFunctionMesh>(
    1.0,                    // 1cm 网格分辨率
    settings::n_batches,    // 最大批次数
    true                    // 自动获取边界
  );
}

// Start new batch for both green function meshes
if (settings::clutch_on) {
  if (simulation::green_function_mesh) {
    simulation::green_function_mesh->start_new_batch(simulation::current_batch);
  }
  if (simulation::fission_green_function_mesh) {
    simulation::fission_green_function_mesh->start_new_batch(simulation::current_batch);
  }
}
```

#### 2.3 修改结束代码

**位置**: `finalize_batch()` 函数中，约第 513-517 行

**原代码**:

```cpp
if (settings::clutch_on && simulation::green_function_mesh) {
  if (simulation::current_batch == settings::n_batches) {
    simulation::green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch);
  }
}
```

**修改为**:

```cpp
// Finalize both green function meshes at the end
if (settings::clutch_on && simulation::current_batch == settings::n_batches) {
  // 输出通量格林函数数据
  if (simulation::green_function_mesh) {
    simulation::green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch);
  }
  
  // 输出裂变源格林函数数据
  if (simulation::fission_green_function_mesh) {
    simulation::fission_green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch);
  }
}
```

---

### 步骤 3: 修改 `src/particle.cpp`

**位置**: `event_advance()` 函数中，第 307-318 行

**原代码**: （保持不变，但添加注释）

```cpp
// 累积通量格林函数贡献
if (simulation::green_function_mesh && material() != MATERIAL_VOID) {
  // 计算通量贡献（权重 × 距离 / 总截面）
  double contribution = 0.0;
  if (macro_xs().total > 0.0) {
    contribution = wgt() * distance / macro_xs().total;
  } else {
    contribution = wgt() * distance; // 如果总截面为0，直接使用权重×距离
  }
  simulation::green_function_mesh->accumulate(
    r(), contribution, source_particle_id());
}
```

**说明**: 这部分代码不需要修改，只是将变量名从 `green_function_mesh` 保持为通量格林函数专用。

---

### 步骤 3: 传递函数统计代码

在 `src/tallies/tally_scoring.cpp` 的 `score_general_ce_nonanalog()` 函数中：

```cpp
case SCORE_GREENFUNCTION:
  // 计算传递函数贡献 (Transfer Function)
  // contribution = 期望裂变中子数 nu_t
  if (p.type() == Type::neutron && p.fission()) {
    // 计算期望裂变中子数（与 physics.cpp 中的 nu_t 相同）
    double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;
    double nu_t = p.wgt() / simulation::keff * weight *
                  p.neutron_xs(p.event_nuclide()).nu_fission /
                  p.neutron_xs(p.event_nuclide()).total;
    
    // 累积到传递函数网格
    if (simulation::transfer_function_mesh && p.source_particle_id() != -1) {
      simulation::transfer_function_mesh->accumulate(
        p.r(), nu_t, p.source_particle_id());
    }
  }
  score = 0.0; // tally 本身不记录数值，只作为开关
  break;
```

**关键点**:

- 使用期望裂变中子数 `nu_t` 作为贡献值
- 与 `physics.cpp` 中裂变矩阵使用相同的计算公式
- 仅在发生裂变事件时统计
- 统计逻辑与 tally 系统集成，便于维护

---

### 步骤 5: 修改 `src/greenfunction_mesh.cpp`

**位置**: `finalize_greenfunction_mesh()` 函数

需要修改输出文件名，以区分两种格林函数：

**原代码**:

```cpp
void GreenFunctionMesh::finalize_greenfunction_mesh(const int batch_id)
{
  // ...
  
  // 创建HDF5文件
  hid_t file_id = file_open("green_function_data.h5", 'w');
  
  // ...
}
```

**问题**: 两个网格都输出到同一个文件会互相覆盖！

**解决方案**: 添加文件名参数或使用不同的默认名称

#### 方案 A: 添加文件名参数（推荐）

修改 `greenfunction_mesh.h`:

```cpp
class GreenFunctionMesh {
public:
  // 修改 finalize 函数签名
  void finalize_greenfunction_mesh(
    const int batch_id,
    const std::string& filename = "green_function_data.h5"
  );
};
```

修改 `greenfunction_mesh.cpp`:

```cpp
void GreenFunctionMesh::finalize_greenfunction_mesh(
  const int batch_id,
  const std::string& filename)
{
  // 保存最后一个batch的数据
  start_new_batch(-1);
  
  if (particle_green_functions_.empty()) {
    return;
  }
  
  // 使用参数指定的文件名
  hid_t file_id = file_open(filename, 'w');
  
  // ... 其余代码不变 ...
}
```

然后在 `simulation.cpp` 中调用时指定不同的文件名:

```cpp
// Finalize both green function meshes with different filenames
if (settings::clutch_on && simulation::current_batch == settings::n_batches) {
  if (simulation::green_function_mesh) {
    simulation::green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch,
      "flux_green_function_data.h5"  // 通量格林函数文件
    );
  }
  
  if (simulation::fission_green_function_mesh) {
    simulation::fission_green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch,
      "fission_green_function_data.h5"  // 裂变源格林函数文件
    );
  }
}
```

---

## 📂 输出文件说明

将生成一个 HDF5 文件：

### `transfer_function_data.h5`

- **内容**: 每个源粒子的传递函数矩阵
- **物理意义**: T(P₀ → r) = 从源点 P₀ 出发的中子在 r 处产生的期望裂变中子数
- **数据大小**: 取决于裂变事件频率（相对较小）
- **统计量**: 期望裂变中子数 nu_t

文件内部结构：

```
transfer_function_data.h5
├── 属性
│   ├── filetype
│   ├── pitch
│   └── n_source_particles
├── 数据集
│   ├── shape [nx, ny, nz]
│   ├── origin [x0, y0, z0]
│   ├── cumulative_green_function [N]
│   └── source_particle_ids [M]
└── 群组
    └── source_particles/
        ├── particle_1 [N]
        ├── particle_10 [N]
        └── ...
```

---

## 🔍 验证方法

### Python 验证脚本

```python
import h5py
import numpy as np

# 读取传递函数数据
with h5py.File('transfer_function_data.h5', 'r') as f:
    # 获取网格参数
    shape = f['shape'][:]
    origin = f['origin'][:]
    pitch = f.attrs['pitch']
    
    # 获取源粒子信息
    source_ids = f['source_particle_ids'][:]
    n_particles = len(source_ids)
    
    print(f"网格尺寸: {shape}")
    print(f"网格原点: {origin}")
    print(f"网格间距: {pitch} cm")
    print(f"源粒子数量: {n_particles}")
    
    # 读取累积传递函数
    cumulative_tf = f['cumulative_green_function'][:]
    total_contribution = np.sum(cumulative_tf)
    
    print(f"\n传递函数总贡献: {total_contribution:.6e}")
    print(f"平均每个网格单元: {total_contribution / np.prod(shape):.6e}")
    
    # 分析单个源粒子
    particle_id = source_ids[0]
    particle_tf = f[f'source_particles/particle_{particle_id}'][:]
    
    print(f"\n源粒子 {particle_id}:")
    print(f"  总贡献: {np.sum(particle_tf):.6e}")
    print(f"  非零网格数: {np.count_nonzero(particle_tf)}")
```

**预期结果**:

- 传递函数仅在发生裂变的区域有非零值
- 每个源粒子的贡献值范围应该在合理的物理范围内
- 总贡献值应该与裂变矩阵的数值量级一致

---

## ⚠️ 注意事项

1. **计算资源**:
   - 内存估算: `N_particles × (nx × ny × nz) × 8 bytes`
   - 建议: 对于大规模计算，考虑使用较粗的网格或减少粒子数

2. **性能影响**:
   - 仅在裂变事件时累积，性能影响小
   - 与裂变矩阵计算同时进行，几乎无额外开销

3. **数据一致性**:
   - 使用与裂变矩阵相同的物理量 `nu_t`
   - 确保 tally 定义正确以启用统计

4. **后处理分析**:
   - 可用于计算伴随通量
   - 可用于重要性函数分析
   - 可用于扰动分析

---

## 📊 应用示例

### 示例 1: 传递函数可视化

```python
import matplotlib.pyplot as plt
import h5py
import numpy as np

with h5py.File('transfer_function_data.h5', 'r') as f:
    # 读取单个源粒子的传递函数
    particle_id = f['source_particle_ids'][0]
    tf_data = f[f'source_particles/particle_{particle_id}'][:]
    
    shape = f['shape'][:]
    tf_3d = tf_data.reshape(shape)
    
    # 绘制中心截面
    fig, ax = plt.subplots(figsize=(10, 8))
    
    center_z = shape[2] // 2
    im = ax.imshow(tf_3d[:, :, center_z].T, cmap='hot')
    ax.set_title(f'传递函数 - 源粒子 {particle_id}')
    ax.set_xlabel('X 网格索引')
    ax.set_ylabel('Y 网格索引')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('期望裂变中子数')
    
    plt.tight_layout()
    plt.savefig('transfer_function_visualization.png', dpi=300)
```

### 示例 2: 源重要性分析

```python
# 计算每个源粒子的总贡献（重要性）
import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File('transfer_function_data.h5', 'r') as f:
    source_ids = f['source_particle_ids'][:]
    
    importances = []
    for sid in source_ids:
        tf = f[f'source_particles/particle_{sid}'][:]
        importances.append(np.sum(tf))
    
    # 绘制重要性分布直方图
    plt.figure(figsize=(12, 6))
    plt.hist(importances, bins=50, edgecolor='black')
    plt.xlabel('总传递函数值（重要性）')
    plt.ylabel('源粒子数量')
    plt.title('源粒子重要性分布')
    plt.grid(True, alpha=0.3)
    plt.savefig('source_importance_distribution.png', dpi=300)
```

---

## 🎯 总结

当前的传递函数实现提供：

✅ **单一物理量** - 期望裂变中子数 `nu_t`，物理意义明确

✅ **Tally 控制** - 通过定义 `<scores>greenfunction</scores>` 启用统计

✅ **代码维护性** - 统计逻辑集中在 tally 系统中，与物理过程分离

✅ **与裂变矩阵一致** - 使用相同的计算公式，确保物理一致性

✅ **高效计算** - 仅在裂变事件时累积，性能影响小

✅ **数据完整性** - 避免混合累积造成的数据污染

---

**下一步**: 请确认是否需要我实际执行这些修改？
