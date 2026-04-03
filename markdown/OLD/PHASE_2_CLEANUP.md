# Phase 2.1 代码清理说明

## 清理日期

2025年11月26日

## 清理内容

Phase 2.2 已完全替代 Phase 2.1 的简化实现，因此清理了所有过时的注释和标记。

### 1. 源文件修改 (`src/beta_effective.cpp`)

#### 简化控制台输出

**修改前:**

```cpp
std::cout << "  Mode: Phase 2 (MATERIAL_DEPENDENT)" << std::endl;
std::cout << "\n  Building cell-material mapping..." << std::endl;
build_cell_material_map(flux);
std::cout << "  Mapped cells: " << cell_to_material_.size() << std::endl;
std::cout << "  Unique fissionable materials: " << unique_materials_.size() << std::endl;
```

**修改后:**

```cpp
std::cout << "  Mode: Phase 2.2 (Geometry-based material mapping)" << std::endl;
build_cell_material_map(flux);
// 详细统计信息已移至 build_cell_material_map() 函数内部
```

#### 清理核数据提取注释

**修改前:**

```cpp
// Phase 2.1 简化: 使用热中子能量的截面 (0.0253 eV)
// TODO Phase 3: 多能群扩展

// 获取裂变截面 (注意: 这需要访问核素的截面数据)
// 简化实现: 假设已知的典型值

// 根据核素名称使用已知数据 (Phase 2.1 简化)
```

**修改后:**

```cpp
// 使用热中子能量的截面 (0.0253 eV)
// TODO Phase 3: 从核数据库提取能量相关截面

// 获取裂变截面
// 当前使用硬编码的热中子典型值，未来将从核数据库动态提取

// 根据核素名称使用已知的热中子数据
```

#### 更新 HDF5 元数据

**修改前:**

```cpp
nuclear_data_source = "Material-dependent (Phase 2.1)";
```

**修改后:**

```cpp
nuclear_data_source = "Material-dependent (Phase 2.2 - Geometry query)";
```

---

### 2. 头文件修改 (`include/openmc/beta_effective.h`)

#### 更新方法注释

**修改前:**

```cpp
//! 构建单元-材料映射 (Phase 2)
//! 通过网格中心点采样确定每个单元的主导材料
void build_cell_material_map(...);

//! 从网格单元索引计算中心位置 (Phase 2.2)
Position grid_index_to_position(int cell_idx) const;
```

**修改后:**

```cpp
//! 构建单元-材料映射 (Phase 2.2)
//! 通过几何查询确定每个网格单元的实际材料
void build_cell_material_map(...);

//! 从网格单元索引计算中心位置
//! (用于几何查询)
Position grid_index_to_position(int cell_idx) const;
```

#### 更新数据成员注释

**修改前:**

```cpp
// Phase 2: 材料相关核数据
std::unordered_map<int, MaterialNuclearData> cell_nuclear_data_;
```

**修改后:**

```cpp
// Phase 2.2: 材料相关核数据 (几何查询)
std::unordered_map<int, MaterialNuclearData> cell_nuclear_data_;
```

---

## 保留的 Phase 1 代码

Phase 1 (FIXED_U235 模式) 的代码**完全保留**，用于：

- 验证和基准测试
- 简单问题的快速计算
- 与 Phase 2.2 结果对比

---

## 编译验证

```bash
cd build
mingw32-make -j24
```

**预期输出:**

- 编译成功，无警告
- `libopenmc.dll` 和 `openmc.exe` 生成

---

## 运行时输出变化

### Phase 2.2 模式输出示例

```
======================================================================
EFFECTIVE DELAYED NEUTRON FRACTION COMPUTATION

[1/4] Reading forward flux from: flux_mesh.h5
  Grid: 12 x 12 x 12
  Pitch: 1.0 cm
  Non-zero cells: 1728

[2/4] Reading adjoint flux from: adjoint_flux.h5
  Grid: 12 x 12 x 12
  Pitch: 1.0 cm
  Non-zero cells: 1728

[3/4] Computing β_eff using MATERIAL-DEPENDENT nuclear data
  Mode: Phase 2.2 (Geometry-based material mapping)
  Cell volume = 1.000 cm³

  Building cell-material mapping using geometry queries...
  Total materials in library: 2
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
  
  Denominator = 1.234567e+05

  Delayed group contributions:
  Group    Numerator        β_i,eff
  --------------------------------------------------
    1    5.123456e+01   4.153204e-04
    2    3.456789e+02   2.801234e-03
    ...

[4/4] Writing results to: beta_eff.h5
  Results written successfully
  β_eff = 0.00267
  Based on 2 fissionable material(s)
======================================================================
```

---

## HDF5 输出变化

### metadata 组

```python
f['metadata'].attrs['beta_eff_mode']           # "MATERIAL_DEPENDENT"
f['metadata'].attrs['nuclear_data_source']     # "Material-dependent (Phase 2.2 - Geometry query)"
```

### diagnostics 组

```python
# Phase 2.2 新增
f['diagnostics/material_cell_counts'][:]       # 每个材料的单元数统计
f['diagnostics'].attrs['total_fissionable_cells']  # 总可裂变单元数
```

---

## 清理效果总结

✅ **移除所有 "Phase 2.1 simplification" 标记**  
✅ **更新所有注释为 Phase 2.2 标准**  
✅ **简化冗余的控制台输出**  
✅ **保持 Phase 1 代码用于验证**  
✅ **HDF5 元数据反映正确的实现版本**  

代码现在更加清晰，准确反映了当前的 Phase 2.2 实现（基于几何查询的材料映射）。
