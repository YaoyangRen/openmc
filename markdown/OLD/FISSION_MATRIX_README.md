# 裂变矩阵功能文档索引

## 📚 文档列表

### 1. [快速参考](./FISSION_MATRIX_QUICK_REFERENCE.md)

**适合**: 快速查询和日常使用

**内容**:

- ✅ 三种运行模式快速切换
- ✅ 编译和验证步骤
- ✅ Python 读取示例
- ✅ 内存使用估算

---

### 2. [运行控制指南](./FISSION_MATRIX_CONTROL_GUIDE.md)

**适合**: 详细配置和优化

**内容**:

- 📖 完整的模式说明
- 📖 批次控制详解
- 📖 性能优化建议
- 📖 常见问题解答
- 📖 网格尺寸推荐

---

### 3. [稀疏矩阵详解](./SPARSE_FISSION_MATRIX.md)

**适合**: 数据分析和后处理

**内容**:

- 🔬 COO 格式规范
- 🔬 Python 完整读取示例
- 🔬 SciPy 稀疏矩阵操作
- 🔬 可视化工具
- 🔬 特征值分析（k-effective）
- 🔬 性能基准测试

---

## 🚀 快速开始

### 第一次使用

1. **阅读**: [快速参考](./FISSION_MATRIX_QUICK_REFERENCE.md) （5分钟）
2. **配置**: 根据需要选择运行模式
3. **编译**: `mingw32-make -j24`
4. **运行**: 执行仿真，获得 `fission_matrix.h5`
5. **分析**: 使用 Python 读取和分析数据

### 常见任务

#### 切换运行模式

→ [快速参考 - 三种模式切换](./FISSION_MATRIX_QUICK_REFERENCE.md#三种模式快速切换)

#### 读取 HDF5 数据

→ [稀疏矩阵详解 - Python 读取](./SPARSE_FISSION_MATRIX.md#python-读取示例)

#### 计算特征值

→ [稀疏矩阵详解 - 特征值分析](./SPARSE_FISSION_MATRIX.md#特征值分析)

#### 优化内存使用

→ [运行控制指南 - 性能考虑](./FISSION_MATRIX_CONTROL_GUIDE.md#性能考虑)

---

## 💡 核心特性

### ✨ 稀疏矩阵存储

**自动优化**: 无需手动配置

- 内存节省: 95%+ （典型稀疏度 99%）
- 格式: COO (Coordinate) 标准格式
- 兼容: SciPy, MATLAB, Eigen

**内存对比**:

| 网格 | 密集 | 稀疏 | 节省 |
|-----|-----|------|------|
| 11³ | 13.5 MB | 0.4 MB | 97% |
| 50³ | 117 GB | 35.7 MB | 99.97% |

### 🎯 灵活的运行控制

**三种模式**:

1. 仅活跃代（推荐）- 精确统计
2. 所有批次 - 完整演化
3. 仅非活跃代 - 收敛调试

### 📊 标准数据格式

**HDF5 输出** (`fission_matrix.h5`):

```
├── 网格信息: origin, shape, pitch
├── 稀疏数据: row_indices, col_indices
├── 原始值: data_raw
└── 归一化值: data_normalized
```

---

## 🔧 技术细节

### 物理原理

裂变矩阵 $F[i][j]$ 定义:

- **行 $i$**: 源粒子出生的空间单元
- **列 $j$**: 裂变事件发生的空间单元
- **值**: 从单元 $i$ 出发导致单元 $j$ 裂变的期望中子数

$$F[i][j] = \sum_{\text{源在}i} \frac{w \cdot \nu \sigma_f}{\sigma_t \cdot k_{\text{eff}}}$$

### 实现细节

- **语言**: C++17
- **并发**: OpenMP + std::mutex
- **存储**: std::unordered_map (COO)
- **输出**: HDF5 with attributes

### 代码位置

- **头文件**: `include/openmc/fission_matrix.h`
- **实现**: `src/fission_matrix.cpp`
- **集成**: `src/simulation.cpp`, `src/particle.cpp`, `src/physics.cpp`

---

## 📖 相关文档

### 格林函数功能

- [GREEN_FUNCTION_GUIDE.md](./GREEN_FUNCTION_GUIDE.md)
- [GREEN_FUNCTION_MODIFICATION_GUIDE.md](./GREEN_FUNCTION_MODIFICATION_GUIDE.md)
- [DUAL_GREEN_FUNCTION_COMPLETION_REPORT.md](./DUAL_GREEN_FUNCTION_COMPLETION_REPORT.md)

### 其他参考

- [SOURCE_TRACKING_README.md](./SOURCE_TRACKING_README.md) - 源粒子追踪

---

## 🆘 获取帮助

### 常见问题

**Q: 稀疏矩阵如何转密集矩阵？**

```python
F_dense = F_sparse.toarray()  # 注意：大矩阵会耗尽内存！
```

**Q: 如何计算 k-effective？**
→ 见 [SPARSE_FISSION_MATRIX.md - 特征值分析](./SPARSE_FISSION_MATRIX.md#特征值分析)

**Q: 输出文件太大？**
→ 稀疏存储已自动压缩，通常文件 < 100 MB

**Q: 稀疏度太低？**
→ 增大网格间距，或检查几何耦合度

### 问题反馈

如有问题，请提供：

1. 网格配置（shape, pitch）
2. 运行模式（活跃代/全部/非活跃代）
3. 统计信息（nnz, sparsity）
4. 错误信息（如有）

---

**文档版本**: v2.0 (稀疏矩阵版本)  
**最后更新**: 2025年10月27日
