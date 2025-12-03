# Phase 2.2 实现说明

## 概述

Phase 2.2 为 β_eff 计算实现了基于真实几何查询的材料映射，替换了 Phase 2.1 中所有网格单元共享单一材料的简化方案。

## 实现日期

2025年11月26日

## 主要变更

### 1. 头文件增强 (`include/openmc/beta_effective.h`)

#### 新增方法

```cpp
// 从一维网格索引计算三维索引
std::array<int, 3> get_grid_indices(int cell_idx) const;

// 从网格单元索引计算中心位置 (Phase 2.2 核心方法)
Position grid_index_to_position(int cell_idx) const;
```

#### 前向声明

```cpp
struct Position;  // 避免包含完整的 position.h
```

---

### 2. 源文件实现 (`src/beta_effective.cpp`)

#### 新增头文件包含

```cpp
#include "openmc/particle_data.h"  // GeometryState
#include "openmc/position.h"       // Position
```

#### 重写 `build_cell_material_map()` - 几何查询核心逻辑

**Phase 2.1 (旧):**

- 使用第一个可裂变材料作为所有单元的默认值
- 简化假设：空间均匀分布

**Phase 2.2 (新):**

```cpp
for (const auto& [cell_idx, flux_val] : flux) {
    // 1. 计算网格单元中心坐标
    Position center_pos = grid_index_to_position(cell_idx);
    
    // 2. 创建几何状态对象
    GeometryState geom;
    geom.r() = center_pos;
    geom.u() = Direction{0.0, 0.0, 1.0};  // 任意方向
    
    // 3. 执行几何查询
    if (!exhaustive_find_cell(geom)) {
        geometry_failed_count++;
        continue;  // 不在几何体内
    }
    
    // 4. 获取材料索引 (处理 lattice 和 distribcell)
    int cell_index = geom.lowest_coord().cell();
    Cell* cell = model::cells[cell_index].get();
    int instance = geom.cell_instance();
    int material_idx = cell->material(instance);
    
    // 5. 检查可裂变性并提取核数据
    if (material_idx != MATERIAL_VOID && material->fissionable()) {
        // 缓存材料核数据
        if (material_data_cache.find(material_id) == material_data_cache.end()) {
            material_data_cache[material_id] = 
                extract_material_nuclear_data(material_id);
        }
        
        // 存储映射
        cell_to_material_[cell_idx] = material_id;
        cell_nuclear_data_[cell_idx] = material_data_cache[material_id];
    }
}
```

**关键特性：**

- ✅ 使用 `exhaustive_find_cell()` 进行真实几何查询
- ✅ 正确处理 lattice 和 distributed cell 的复杂几何
- ✅ 材料核数据缓存机制（避免重复提取）
- ✅ 统计 void、非裂变、几何查询失败的单元

---

#### 新增辅助方法实现

**`get_grid_indices(int cell_idx)`**

```cpp
// 从一维索引反推三维网格索引 (row-major: z,y,x)
int iz = cell_idx / (nx * ny);
int remainder = cell_idx % (nx * ny);
int iy = remainder / nx;
int ix = remainder % nx;
return {ix, iy, iz};
```

**`grid_index_to_position(int cell_idx)`**

```cpp
auto [ix, iy, iz] = get_grid_indices(cell_idx);

// 计算网格单元中心的全局坐标
double x = grid_lower_left_[0] + (ix + 0.5) * grid_pitch_;
double y = grid_lower_left_[1] + (iy + 0.5) * grid_pitch_;
double z = grid_lower_left_[2] + (iz + 0.5) * grid_pitch_;

return Position{x, y, z};
```

---

### 3. HDF5 诊断数据增强

**新增输出数据集 (diagnostics 组):**

```python
# Phase 2.2 材料分布统计
material_ids              # [n_materials] 材料ID列表
material_nu_total         # [n_materials] 各材料的ν_total
material_sigma_f          # [n_materials] 各材料的宏观裂变截面
material_cell_counts      # [n_materials] 各材料的单元数 (NEW!)

# 属性
n_unique_materials        # 唯一可裂变材料数
total_fissionable_cells   # 总可裂变单元数 (NEW!)
```

**Python 读取示例：**

```python
import h5py

with h5py.File('beta_eff.h5', 'r') as f:
    mat_ids = f['diagnostics/material_ids'][:]
    mat_counts = f['diagnostics/material_cell_counts'][:]
    
    for mid, count in zip(mat_ids, mat_counts):
        print(f"Material {mid}: {count} cells")
```

---

### 4. 控制台输出增强

**Phase 2.2 新增诊断信息：**

```
Building cell-material mapping using geometry queries...
  Total materials in library: 3
  Querying geometry for 1728 mesh cells...
  ------------------------------------------------------------
  Geometry query results:
    Fissionable cells: 1500
    Non-fissionable cells: 100
    Void cells: 50
    Geometry query failed: 78
    Unique fissionable materials: 2

  Material distribution:
    Material 1 (fuel): 800 cells
      ν_total = 2.4300, Σ_f = 5.840e-01 cm⁻¹
    Material 2 (blanket): 700 cells
      ν_total = 2.8800, Σ_f = 7.474e-01 cm⁻¹
  ------------------------------------------------------------
```

---

## 技术细节

### 几何查询 API

**`exhaustive_find_cell(GeometryState& geom)`**

- **功能**: 在整个几何树中定位粒子位置
- **返回**: `true` 成功找到，`false` 不在几何体内
- **副作用**: 填充 `geom` 的坐标层级、cell、材料索引等信息
- **线程安全**: 安全（每个线程使用独立的 `GeometryState` 对象）

**`Cell::material(int instance)`**

- **功能**: 获取指定实例的材料索引
- **参数**: `instance` - cell 实例编号（处理 lattice 和 distribcell）
- **返回**: 材料索引（非 ID！）

### 坐标系统

**网格坐标系:**

- 原点: `grid_lower_left_` (从 flux_mesh.h5 读取)
- 间距: `grid_pitch_` (均匀网格)
- 索引顺序: row-major (z, y, x)

**一维索引到三维映射:**

```
cell_idx = iz * (nx * ny) + iy * nx + ix
```

**网格中心坐标:**

```
center = lower_left + (index + 0.5) * pitch
```

---

## 验证测试

### 测试场景 1: 均匀燃料堆芯

```xml
<material id="1" name="fuel">
  <nuclide name="U235" ao="0.02" />
</material>
```

**预期结果:**

- 所有单元映射到材料 1
- `n_unique_materials = 1`
- `material_cell_counts = [total_flux_cells]`

---

### 测试场景 2: 燃料 + 反射层

```xml
<material id="1" name="fuel">
  <nuclide name="U235" ao="0.02" />
</material>
<material id="2" name="reflector">
  <nuclide name="H1" ao="0.06" />
</material>
```

**预期结果:**

- 燃料区单元 → 材料 1
- 反射层单元 → 跳过（非裂变）
- `n_unique_materials = 1`
- 输出显示 "Non-fissionable cells: N"

---

### 测试场景 3: Pin-cell Lattice

```xml
<cell id="1" material="1" region="..." />  <!-- 燃料棒 -->
<cell id="2" material="2" region="..." />  <!-- 慢化剂 -->
<lattice id="10" type="rect">
  <pitch>1.26 1.26</pitch>
  <universes>1 1 1 ...</universes>
</lattice>
```

**预期结果:**

- 燃料棒位置单元 → 材料 1
- 慢化剂位置单元 → 跳过（通常非裂变）
- 正确处理 lattice instance 编号

---

## 性能考虑

### 计算复杂度

- **Phase 2.1**: O(1) - 单次材料提取
- **Phase 2.2**: O(n_cells × geometry_depth)
  - `n_cells`: 通量网格单元数（~1000-10000）
  - `geometry_depth`: 几何树深度（~3-10）

### 优化策略

1. **材料核数据缓存**:
   - 每个唯一材料只提取一次核数据
   - 缓存在 `material_data_cache` map 中

2. **早期跳过**:
   - 几何查询失败 → 立即 continue
   - MATERIAL_VOID → 立即 continue
   - 非裂变材料 → 立即 continue

3. **未来优化方向**:
   - 并行化几何查询（OpenMP）
   - 空间局部性优化（按 z-order 排序）

---

## 局限性与未来工作

### 当前局限性 (Phase 2.2)

1. **单点采样**
   - 使用网格单元中心点
   - 对于跨越材料边界的单元可能不准确

2. **体积权重缺失**
   - 未考虑部分体积占用
   - 假设中心材料代表整个单元

3. **硬编码热中子截面**
   - σ_f, ν 使用固定值（0.0253 eV）
   - 未从核数据库动态提取

---

### Phase 2.3 计划：多点采样

**目标**: 提高跨材料边界单元的精度

**方法**:

```cpp
// 8点分层采样 (立方体顶点)
for (int s = 0; s < 8; ++s) {
    Position sample = lower + 
        Position((s&1)*pitch, ((s>>1)&1)*pitch, ((s>>2)&1)*pitch);
    
    // 几何查询
    if (exhaustive_find_cell(...)) {
        material_hits[mat_idx]++;
    }
}

// 加权累加
double weight = material_hits[mat_id] / 8.0;
```

---

### Phase 3 计划：真实截面库

**目标**: 从 OpenMC 核数据库提取能量相关截面

**API 调用**:

```cpp
// 从 Nuclide 对象获取裂变反应
const auto& fission_rxn = nuc->reactions_[FISSION];

// 在热中子能量点插值
double E_thermal = 0.0253e-6;  // MeV
double sigma_f = fission_rxn.xs_->operator()(E_thermal);
```

---

## 编译与测试

### 编译

```bash
cd build
mingw32-make -j24
```

### 运行测试

```bash
# 使用 MATERIAL_DEPENDENT 模式
./bin/openmc

# 检查 beta_eff.h5 输出
python -c "
import h5py
with h5py.File('beta_eff.h5', 'r') as f:
    print('Mode:', f['metadata'].attrs['beta_eff_mode'])
    print('Unique materials:', f['diagnostics'].attrs['n_unique_materials'])
    print('Material IDs:', f['diagnostics/material_ids'][:])
    print('Cell counts:', f['diagnostics/material_cell_counts'][:])
"
```

---

## 代码示例

### 完整使用流程

```cpp
// 1. 创建 BetaEffective 对象 (Phase 2.2 模式)
BetaEffective beta_calc; // 默认启用材料相关模式

// 2. 从文件计算 β_eff
beta_calc.compute_from_files(
    "flux_mesh.h5", 
    "adjoint_flux.h5", 
    "beta_eff.h5"
);

// 3. 读取结果
double beta_total = beta_calc.get_beta_total();
auto beta_i = beta_calc.get_all_beta_i();

std::cout << "β_eff = " << beta_total << std::endl;
for (int i = 0; i < 6; ++i) {
    std::cout << "  β_" << (i+1) << " = " << beta_i[i] << std::endl;
}
```

---

## 故障排查

### 常见问题

**Q1: "Geometry query failed: N" 数量很大**

- **原因**: 网格超出几何边界
- **解决**: 检查 `grid_lower_left` 和网格尺寸是否正确

**Q2: 所有单元都是 void**

- **原因**: 坐标系不匹配
- **检查**: `flux_mesh.h5` 中的 `grid_lower_left` 是否与几何一致

**Q3: 编译错误 "Position not declared"**

- **原因**: 缺少前向声明
- **解决**: 确保 `beta_effective.h` 中有 `struct Position;`

---

## 总结

Phase 2.2 实现了**真实的几何查询驱动材料映射**，显著提升了 β_eff 计算的物理保真度：

✅ **空间变化材料支持**: 每个网格单元查询实际材料分布  
✅ **复杂几何兼容**: 正确处理 lattice、distribcell、nested universes  
✅ **诊断数据完善**: 输出材料分布统计，便于验证  
✅ **代码健壮性**: 完善的错误检查和边界情况处理  

**下一步**: Phase 2.3 将实现多点采样，提高跨材料边界单元的精度。
