# 双格林函数网格修改指南

本文档展示如何同时保留**通量格林函数**和**裂变源格林函数**两种贡献，并将它们分开存储。

---

## 📊 两种格林函数的物理意义

### 1. 通量格林函数 (Flux Green's Function)
- **物理意义**: 单位点源在空间位置 r 处产生的中子通量
- **计算公式**: `Φ(r | r') = ∫ w(s) × ds / Σ_t`
- **累积时机**: 粒子每次移动时（`event_advance()`）
- **应用**: 源重要性分析、探测器响应计算

### 2. 裂变源格林函数 (Fission Source Green's Function)
- **物理意义**: 单位点源产生的裂变中子在空间的分布
- **计算公式**: `S_f(r | r') = w × ν × Σ_f / Σ_t`
- **累积时机**: 发生裂变碰撞时（`event_collide()` → `score_collision_tally()`）
- **应用**: 裂变链反应分析、临界安全评估

---

## 🔧 修改步骤

### 步骤 1: 修改 `include/openmc/simulation.h`

```cpp
namespace simulation {

// 现有变量...
extern "C" int current_batch;
extern "C" int current_gen;

// 添加两个格林函数网格
extern std::unique_ptr<GreenFunctionMesh> 
  green_function_mesh;          // 通量格林函数 (Flux)
  
extern std::unique_ptr<GreenFunctionMesh> 
  fission_green_function_mesh;  // 裂变源格林函数 (Fission Source)

} // namespace simulation
```

---

### 步骤 2: 修改 `src/simulation.cpp`

#### 2.1 添加全局变量声明

**位置**: 约第 326 行，在 `green_function_mesh` 声明后添加

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

### 步骤 4: 修改 `src/tallies/tally_scoring.cpp`

**位置**: `score_collision_tally()` 函数中，约第 1003-1018 行

**原代码**:
```cpp
case SCORE_GREENFUNCTION:
  if (settings::clutch_on) {
    if (p.type() == Type::neutron && (p.fission())) {
      double contribution = 0.0;
      if (p.neutron_xs(p.event_nuclide()).total > 0) {
        contribution =
          p.wgt_last() * p.neutron_xs(p.event_nuclide()).nu_fission *
          p.neutron_xs(p.event_nuclide()).fission /
          p.neutron_xs(p.event_nuclide()).total; // TODO 是否需要除sigma_t
      }
      if (simulation::green_function_mesh) {
        simulation::green_function_mesh->accumulate(
          p.r(), contribution, p.source_particle_id());
      }
    }
  }
  break;
```

**修改为**:
```cpp
case SCORE_GREENFUNCTION:
  if (settings::clutch_on) {
    if (p.type() == Type::neutron && (p.fission())) {
      // 计算裂变源贡献
      double contribution = 0.0;
      if (p.neutron_xs(p.event_nuclide()).total > 0) {
        contribution =
          p.wgt_last() * p.neutron_xs(p.event_nuclide()).nu_fission *
          p.neutron_xs(p.event_nuclide()).fission /
          p.neutron_xs(p.event_nuclide()).total;
      }
      
      // 累积到裂变源格林函数网格（注意这里改为 fission_green_function_mesh）
      if (simulation::fission_green_function_mesh) {
        simulation::fission_green_function_mesh->accumulate(
          p.r(), contribution, p.source_particle_id());
      }
    }
  }
  break;
```

**关键变化**: 将 `simulation::green_function_mesh` 改为 `simulation::fission_green_function_mesh`

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

修改后将生成两个 HDF5 文件：

### 1. `flux_green_function_data.h5`
- **内容**: 每个源粒子的通量格林函数矩阵
- **物理意义**: G_flux(r, r') = 源位置 r' 的单位点源在 r 处产生的通量
- **数据大小**: 较大（每次粒子移动都累积）

### 2. `fission_green_function_data.h5`
- **内容**: 每个源粒子的裂变源格林函数矩阵
- **物理意义**: G_fission(r, r') = 源位置 r' 产生的裂变中子在 r 处的分布
- **数据大小**: 较小（仅裂变事件累积）

两个文件的内部结构相同：
```
green_function_data.h5
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

# 读取两个文件
with h5py.File('flux_green_function_data.h5', 'r') as f_flux, \
     h5py.File('fission_green_function_data.h5', 'r') as f_fission:
    
    # 验证网格参数一致
    assert np.array_equal(f_flux['shape'][:], f_fission['shape'][:])
    assert np.array_equal(f_flux['origin'][:], f_fission['origin'][:])
    
    # 验证源粒子ID一致
    flux_ids = set(f_flux['source_particle_ids'][:])
    fission_ids = set(f_fission['source_particle_ids'][:])
    
    print(f"通量格林函数粒子数: {len(flux_ids)}")
    print(f"裂变源格林函数粒子数: {len(fission_ids)}")
    print(f"共同粒子数: {len(flux_ids & fission_ids)}")
    
    # 比较累积数据
    flux_total = np.sum(f_flux['cumulative_green_function'][:])
    fission_total = np.sum(f_fission['cumulative_green_function'][:])
    
    print(f"\n通量格林函数总和: {flux_total:.6e}")
    print(f"裂变源格林函数总和: {fission_total:.6e}")
    print(f"比值 (Flux/Fission): {flux_total/fission_total:.2f}")
```

**预期结果**:
- 通量格林函数总和 >> 裂变源格林函数总和（因为累积频率不同）
- 裂变源格林函数仅在发生裂变的区域有非零值
- 通量格林函数在整个几何体中都有分布

---

## ⚠️ 注意事项

1. **内存使用**: 两个网格会占用双倍内存
   - 估算: `2 × N_particles × (nx × ny × nz) × 8 bytes`
   - 建议: 对于大规模计算，考虑使用较粗的网格或减少粒子数

2. **性能影响**: 
   - 通量格林函数累积频率高（每次移动）
   - 裂变源格林函数累积频率低（仅裂变事件）
   - 总体性能影响约 5-10%

3. **数据一致性**:
   - 两个网格的空间参数（shape, origin, pitch）应该相同
   - 建议在初始化时使用相同的参数

4. **后处理分析**:
   - 可以计算两者的比值来分析裂变增殖效应
   - 通量/裂变源比值可以反映系统的增殖能力

---

## 📊 应用示例

### 示例 1: 源重要性对比

```python
# 对比同一个源粒子的两种贡献
import matplotlib.pyplot as plt

with h5py.File('flux_green_function_data.h5', 'r') as f_flux, \
     h5py.File('fission_green_function_data.h5', 'r') as f_fission:
    
    particle_id = 1
    
    flux_data = f_flux[f'source_particles/particle_{particle_id}'][:]
    fission_data = f_fission[f'source_particles/particle_{particle_id}'][:]
    
    shape = f_flux['shape'][:]
    flux_3d = flux_data.reshape(shape)
    fission_3d = fission_data.reshape(shape)
    
    # 绘制中心截面对比
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    center_z = shape[2] // 2
    
    im1 = ax1.imshow(flux_3d[:, :, center_z].T, cmap='viridis')
    ax1.set_title(f'通量格林函数 - 源粒子 {particle_id}')
    plt.colorbar(im1, ax=ax1)
    
    im2 = ax2.imshow(fission_3d[:, :, center_z].T, cmap='hot')
    ax2.set_title(f'裂变源格林函数 - 源粒子 {particle_id}')
    plt.colorbar(im2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('green_function_comparison.png', dpi=300)
```

### 示例 2: 裂变增殖因子分析

```python
# 计算空间相关的裂变增殖因子
multiplication_factor = flux_3d / (fission_3d + 1e-10)  # 避免除零

plt.figure(figsize=(10, 8))
plt.imshow(multiplication_factor[:, :, center_z].T, cmap='coolwarm')
plt.colorbar(label='通量/裂变源比值')
plt.title('空间相关裂变增殖因子')
plt.savefig('multiplication_factor.png', dpi=300)
```

---

## 🎯 总结

通过以上修改，您将获得：

✅ **通量格林函数** (`flux_green_function_data.h5`)
  - 反映每个源粒子对空间通量的贡献
  - 用于探测器响应、源重要性分析

✅ **裂变源格林函数** (`fission_green_function_data.h5`)
  - 反映每个源粒子产生的裂变链
  - 用于临界安全、裂变增殖分析

✅ **物理意义清晰** - 两种不同的物理量分开存储，便于分析

✅ **数据完整性** - 避免混合累积造成的数据污染

---

**下一步**: 请确认是否需要我实际执行这些修改？
