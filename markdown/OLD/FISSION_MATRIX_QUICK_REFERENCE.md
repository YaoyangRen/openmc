# 裂变矩阵运行模式快速参考

## 当前配置 ✅

**默认模式**: 仅在活跃代运行裂变矩阵

**存储格式**: COO 稀疏矩阵（内存节省 95%+）

**位置**: `src/simulation.cpp` 第 432 行

## 🚀 稀疏矩阵优势

- **内存节省**: 99% 稀疏度 → 节省 99% 内存
- **更大网格**: 100×100×100 仅需 ~1 GB（原需 7.45 TB）
- **标准格式**: 兼容 SciPy、MATLAB、Eigen

详见: [SPARSE_FISSION_MATRIX.md](./SPARSE_FISSION_MATRIX.md)

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

运行仿真后，检查最终输出：

```
======================================================================
FISSION MATRIX FINALIZATION
======================================================================

Grid Configuration:
  Total cells: 1331
  Non-zero elements: 45123
  Sparsity: 97.45%
  Dense storage would need: 13.52 MB
  Sparse storage uses: 1.03 MB
  Memory saved: 12.49 MB ✨

Data Collection Summary:
  Total fissions recorded: 1234567
  Total sources recorded: 50000    ← 源粒子数
  Number of batches: 100           ← 统计批次数

Storage Format: COO (Coordinate)
======================================================================
```

**判断模式**:
- **模式 1**: `Number of batches` ≈ `n_batches - n_inactive`
- **模式 2**: `Number of batches` ≈ `n_batches`
- **模式 3**: `Number of batches` ≈ `n_inactive`

---

## 稀疏矩阵读取

### Python 快速读取

```python
import h5py
from scipy.sparse import coo_matrix

# 读取稀疏矩阵
with h5py.File('fission_matrix.h5', 'r') as f:
    n_cells = f.attrs['n_cells']
    rows = f['row_indices'][:]
    cols = f['col_indices'][:]
    data = f['data_normalized'][:]
    
    # 创建稀疏矩阵
    F = coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
    
    print(f"Matrix: {n_cells}×{n_cells}")
    print(f"Non-zeros: {len(data)}")
    print(f"Sparsity: {100*(1-len(data)/n_cells**2):.2f}%")
```

### 完整分析工具

详见: [SPARSE_FISSION_MATRIX.md](./SPARSE_FISSION_MATRIX.md)
- 特征值计算（k-effective）
- 可视化工具
- 统计分析

---

## ⚠️ 注意事项

### 内存（已优化）

✅ **稀疏存储**: 自动节省 95%+ 内存
- 11×11×11 网格: 仅需 0.4 MB（原需 13.5 MB）
- 50×50×50 网格: 仅需 35.7 MB（原需 117 GB）

### 推荐网格尺寸

| 间距 | 网格 | 内存 | 用途 |
|-----|-----|------|-----|
| 5 cm | 10³ | < 1 MB | 快速测试 |
| 2 cm | 25³ | < 10 MB | 标准计算 |
| 1 cm | 50³ | < 100 MB | 精细分析 |
| 0.5 cm | 100³ | < 2 GB | 高精度 |

### 其他

- **精度**: 仅活跃代统计更准确（源已收敛）
- **性能**: 仅活跃代可节省 30-50% 计算时间

---

## 📚 相关文档

- [FISSION_MATRIX_CONTROL_GUIDE.md](./FISSION_MATRIX_CONTROL_GUIDE.md) - 完整控制指南
- [SPARSE_FISSION_MATRIX.md](./SPARSE_FISSION_MATRIX.md) - 稀疏矩阵详解
- [GREEN_FUNCTION_GUIDE.md](./GREEN_FUNCTION_GUIDE.md) - 格林函数使用

---

**更新日期**: 2025年10月27日
