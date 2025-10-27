# 裂变矩阵运行模式快速参考

## 当前配置 ✅

**默认模式**: 仅在活跃代运行裂变矩阵

**位置**: `src/simulation.cpp` 第 432 行

## 三种模式快速切换

### 📌 模式 1: 仅活跃代（推荐）

```cpp
// 当前激活 ✓
if (settings::clutch_on && simulation::current_batch > settings::n_inactive) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

✅ **用途**: 正常计算，统计收敛后的裂变分布  
📊 **统计**: 批次 [n_inactive+1, n_batches]

---

### 📌 模式 2: 所有批次

```cpp
// 取消下面注释，注释掉模式 1
if (settings::clutch_on) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

✅ **用途**: 研究源收敛过程  
📊 **统计**: 批次 [1, n_batches]

---

### 📌 模式 3: 仅非活跃代

```cpp
// 取消下面注释，注释掉模式 1
if (settings::clutch_on && simulation::current_batch <= settings::n_inactive) {
  if (simulation::fission_matrix) {
    simulation::fission_matrix->start_new_batch(simulation::current_batch);
  }
}
```

✅ **用途**: 调试源收敛问题  
📊 **统计**: 批次 [1, n_inactive]

---

## 示例场景

### 场景: 50 个非活跃代 + 100 个活跃代

```xml
<!-- settings.xml -->
<batches>150</batches>
<inactive>50</inactive>
```

| 模式 | 运行批次 | 统计批次数 | 说明 |
|------|---------|-----------|------|
| 模式 1 | 51-150 | 100 | ✅ 推荐：仅统计收敛后 |
| 模式 2 | 1-150 | 150 | 📊 完整演化过程 |
| 模式 3 | 1-50 | 50 | 🐛 调试用 |

---

## 修改步骤

1. **编辑文件**: `d:\OpenMC\openmc\src\simulation.cpp`
2. **找到位置**: 第 428-448 行（搜索 "Start new batch for fission matrix"）
3. **切换代码**: 注释当前模式，取消注释目标模式
4. **重新编译**: 
   ```powershell
   cd d:\OpenMC\openmc\build
   mingw32-make -j24
   ```

---

## 验证当前模式

运行仿真后，检查输出：

```
Fission matrix written to: fission_matrix.h5
Statistics:
  Total fissions recorded: 1234567
  Total sources recorded: 50000    ← 源粒子数
  ...
```

- **模式 1**: `total_sources` ≈ 活跃代粒子总数
- **模式 2**: `total_sources` ≈ 所有批次粒子总数  
- **模式 3**: `total_sources` ≈ 非活跃代粒子总数

---

## ⚠️ 注意事项

1. **内存**: 网格越细，内存需求越大（∝ n_cells²）
2. **精度**: 仅活跃代统计更准确（源已收敛）
3. **性能**: 仅活跃代可节省 30-50% 计算时间

---

详细说明请参阅: `FISSION_MATRIX_CONTROL_GUIDE.md`
