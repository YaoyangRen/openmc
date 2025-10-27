# 稀疏裂变矩阵存储格式说明

## 概述

裂变矩阵采用 **COO (Coordinate) 格式**的稀疏存储，可以大幅减少内存使用。对于大多数实际问题，裂变矩阵的稀疏度通常超过 95%，意味着内存节省超过 95%。

## 内存节省示例

### 密集存储 vs 稀疏存储

| 网格尺寸 | 总单元数 | 密集矩阵内存 | 稀疏度 | 非零元素 | 稀疏存储内存 | 节省 |
|---------|---------|-------------|--------|---------|-------------|------|
| 10×10×10 | 1,000 | 7.6 MB | 99% | 10,000 | 0.23 MB | **97%** |
| 20×20×20 | 8,000 | 488 MB | 99.5% | 320,000 | 7.3 MB | **98.5%** |
| 50×50×50 | 125,000 | 117 GB | 99.9% | 1,562,500 | 35.7 MB | **99.97%** |
| 100×100×100 | 1,000,000 | 7.45 TB | 99.95% | 50,000,000 | 1.1 GB | **99.985%** |

## HDF5 文件格式

### 文件属性

```python
filetype = "fission_matrix_sparse"
version = "2.0"
storage_format = "COO"  # Coordinate format
pitch = 1.0             # 网格间距 (cm)
n_realizations = 50     # 统计批次数
total_fissions = 123456 # 总裂变事件数
total_sources = 50000   # 总源粒子数
n_cells = 1331          # 总单元数 (nx * ny * nz)
nnz = 45123             # 非零元素数量
```

### 数据集结构

#### 网格信息
- **`origin`** (3,): `[x_min, y_min, z_min]` - 网格原点
- **`shape`** (3,): `[nx, ny, nz]` - 各维度单元数
- **`source_counts`** (n_cells,): 每个单元的源粒子计数

#### 稀疏矩阵数据 (COO格式)
- **`row_indices`** (nnz,): 行索引数组（源单元）
- **`col_indices`** (nnz,): 列索引数组（裂变单元）
- **`data_raw`** (nnz,): 原始裂变中子数
- **`data_normalized`** (nnz,): 归一化后的值

## Python 读取示例

### 方法 1: 使用 h5py 读取

```python
import h5py
import numpy as np

# 打开文件
with h5py.File('fission_matrix.h5', 'r') as f:
    # 读取属性
    n_cells = f.attrs['n_cells']
    nnz = f.attrs['nnz']
    storage_format = f.attrs['storage_format']
    
    print(f"Matrix size: {n_cells} × {n_cells}")
    print(f"Non-zero elements: {nnz}")
    print(f"Sparsity: {100 * (1 - nnz/(n_cells**2)):.2f}%")
    
    # 读取网格信息
    origin = f['origin'][:]
    shape = f['shape'][:]
    
    # 读取稀疏数据
    rows = f['row_indices'][:]
    cols = f['col_indices'][:]
    data_norm = f['data_normalized'][:]
    
    print(f"\nGrid shape: {shape}")
    print(f"Grid origin: {origin}")
```

### 方法 2: 转换为 SciPy 稀疏矩阵

```python
from scipy.sparse import coo_matrix
import h5py

def load_sparse_fission_matrix(filename):
    """加载稀疏裂变矩阵为 SciPy COO 格式"""
    with h5py.File(filename, 'r') as f:
        n_cells = f.attrs['n_cells']
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        
        # 创建稀疏矩阵
        matrix = coo_matrix((data, (rows, cols)), 
                           shape=(n_cells, n_cells))
        
        return matrix

# 使用
F = load_sparse_fission_matrix('fission_matrix.h5')

# 转换为其他稀疏格式
F_csr = F.tocsr()  # Compressed Sparse Row (高效行操作)
F_csc = F.tocsc()  # Compressed Sparse Column (高效列操作)

# 矩阵运算
eigenvalues, eigenvectors = sparse.linalg.eigs(F_csr, k=1)
k_eff = eigenvalues[0].real
```

### 方法 3: 转换为密集矩阵（小矩阵）

```python
import h5py
import numpy as np

def load_as_dense(filename):
    """将稀疏矩阵转换为密集格式（仅适用于小矩阵！）"""
    with h5py.File(filename, 'r') as f:
        n_cells = f.attrs['n_cells']
        
        # 警告：大矩阵会耗尽内存！
        if n_cells > 10000:
            raise ValueError(f"Matrix too large ({n_cells}×{n_cells}). Use sparse format!")
        
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        
        # 初始化密集矩阵
        matrix = np.zeros((n_cells, n_cells))
        matrix[rows, cols] = data
        
        return matrix

# 使用（仅小矩阵）
F_dense = load_as_dense('fission_matrix.h5')
```

## 可视化示例

### 稀疏模式可视化

```python
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import h5py

def visualize_sparsity_pattern(filename):
    """可视化稀疏模式"""
    with h5py.File(filename, 'r') as f:
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        n_cells = f.attrs['n_cells']
    
    plt.figure(figsize=(10, 10))
    plt.spy(coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells)), 
            markersize=1, color='blue')
    plt.title('Fission Matrix Sparsity Pattern')
    plt.xlabel('Fission Cell Index')
    plt.ylabel('Source Cell Index')
    plt.tight_layout()
    plt.savefig('sparsity_pattern.png', dpi=300)
    plt.show()

visualize_sparsity_pattern('fission_matrix.h5')
```

### 热图可视化（小矩阵）

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
import h5py

def plot_fission_matrix_heatmap(filename, max_size=100):
    """绘制裂变矩阵热图"""
    with h5py.File(filename, 'r') as f:
        n_cells = f.attrs['n_cells']
        
        if n_cells > max_size:
            print(f"Matrix too large ({n_cells}), showing {max_size}×{max_size} subset")
            n_cells = max_size
        
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        
        # 筛选子集
        mask = (rows < max_size) & (cols < max_size)
        rows = rows[mask]
        cols = cols[mask]
        data = data[mask]
        
        # 转换为密集矩阵
        matrix = np.zeros((n_cells, n_cells))
        matrix[rows, cols] = data
    
    plt.figure(figsize=(12, 10))
    plt.imshow(matrix, cmap='hot', interpolation='nearest', 
               norm=plt.matplotlib.colors.LogNorm(vmin=1e-6, vmax=matrix.max()))
    plt.colorbar(label='Normalized Fission Rate')
    plt.title('Fission Matrix (log scale)')
    plt.xlabel('Fission Cell Index')
    plt.ylabel('Source Cell Index')
    plt.tight_layout()
    plt.savefig('fission_matrix_heatmap.png', dpi=300)
    plt.show()

plot_fission_matrix_heatmap('fission_matrix.h5')
```

## 特征值分析

### 计算主导特征值（k-effective）

```python
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs
import h5py

def calculate_keff_from_fission_matrix(filename):
    """从裂变矩阵计算 k-effective"""
    with h5py.File(filename, 'r') as f:
        n_cells = f.attrs['n_cells']
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        
        # 创建稀疏矩阵
        F = coo_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
        F_csr = F.tocsr()
    
    # 计算主导特征值
    eigenvalues, eigenvectors = eigs(F_csr, k=1, which='LM')
    k_eff = eigenvalues[0].real
    
    print(f"Dominant eigenvalue (k-eff): {k_eff:.6f}")
    
    # 归一化特征向量（功率分布）
    flux = np.abs(eigenvectors[:, 0])
    flux /= flux.sum()
    
    return k_eff, flux

k_eff, power_dist = calculate_keff_from_fission_matrix('fission_matrix.h5')
```

## 统计分析

### 分析裂变源分布

```python
import h5py
import numpy as np

def analyze_fission_distribution(filename):
    """分析裂变源空间分布"""
    with h5py.File(filename, 'r') as f:
        rows = f['row_indices'][:]
        cols = f['col_indices'][:]
        data = f['data_normalized'][:]
        source_counts = f['source_counts'][:]
        shape = f['shape'][:]
        n_cells = f.attrs['n_cells']
    
    # 计算每个单元的总裂变贡献
    fission_by_cell = np.zeros(n_cells)
    for i, val in zip(cols, data):
        fission_by_cell[i] += val
    
    # 统计
    print("Fission Distribution Analysis:")
    print(f"  Total source particles: {source_counts.sum():.0f}")
    print(f"  Cells with sources: {np.count_nonzero(source_counts)}")
    print(f"  Cells with fissions: {np.count_nonzero(fission_by_cell)}")
    print(f"  Max fission rate: {fission_by_cell.max():.6e}")
    print(f"  Mean fission rate: {fission_by_cell[fission_by_cell>0].mean():.6e}")
    
    # 找到热点
    hot_cells = np.argsort(fission_by_cell)[-10:][::-1]
    print("\nTop 10 fission cells:")
    for rank, cell_idx in enumerate(hot_cells, 1):
        iz = cell_idx // (shape[0] * shape[1])
        iy = (cell_idx % (shape[0] * shape[1])) // shape[0]
        ix = cell_idx % shape[0]
        print(f"  {rank}. Cell [{ix},{iy},{iz}]: {fission_by_cell[cell_idx]:.6e}")

analyze_fission_distribution('fission_matrix.h5')
```

## 内存优化建议

### 1. 网格分辨率选择

```python
# 估算内存需求
def estimate_memory(nx, ny, nz, sparsity=0.99):
    """估算稀疏矩阵内存需求"""
    n_cells = nx * ny * nz
    nnz = int(n_cells * n_cells * (1 - sparsity))
    
    # 每个非零元素需要: 2×int + 1×double ≈ 24 bytes
    sparse_mem_mb = nnz * 24 / (1024**2)
    dense_mem_mb = n_cells * n_cells * 8 / (1024**2)
    
    print(f"Grid: {nx}×{ny}×{nz} = {n_cells} cells")
    print(f"Sparsity: {sparsity*100}%")
    print(f"Non-zero elements: {nnz:,}")
    print(f"Sparse storage: {sparse_mem_mb:.2f} MB")
    print(f"Dense storage: {dense_mem_mb:.2f} MB")
    print(f"Memory saved: {100*(1-sparse_mem_mb/dense_mem_mb):.2f}%")

# 测试不同网格
estimate_memory(50, 50, 50, sparsity=0.995)
```

### 2. 推荐网格尺寸

| 计算类型 | 推荐网格间距 | 典型尺寸 | 内存需求 |
|---------|------------|---------|----------|
| 快速测试 | 5 cm | 10×10×10 | < 1 MB |
| 标准计算 | 2 cm | 25×25×25 | < 10 MB |
| 精细分析 | 1 cm | 50×50×50 | < 100 MB |
| 高精度 | 0.5 cm | 100×100×100 | < 2 GB |

## 性能基准

稀疏矩阵操作性能（相比密集矩阵）：

| 操作 | 密集矩阵 | 稀疏矩阵 (99% 稀疏) | 加速比 |
|------|---------|-------------------|--------|
| 矩阵构建 | O(n²) | O(nnz) | ~100× |
| 内存使用 | O(n²) | O(nnz) | ~100× |
| 矩阵-向量乘 | O(n²) | O(nnz) | ~100× |
| 特征值计算 | O(n³) | O(n×nnz) | ~10-50× |

## 常见问题

### Q1: 如何确定合适的网格间距？

根据几何特征尺寸选择：
- `pitch ≈ 特征长度 / 10` (粗网格)
- `pitch ≈ 特征长度 / 20` (标准)
- `pitch ≈ 特征长度 / 50` (精细)

### Q2: 稀疏度低于 90% 怎么办？

- 增大网格间距
- 检查几何是否高度耦合
- 考虑使用更高效的稀疏格式（CSR/CSC）

### Q3: 如何与其他工具集成？

稀疏矩阵 COO 格式是标准格式，支持：
- SciPy (Python)
- MATLAB (sparse)
- Eigen (C++)
- PETSc (并行计算)

## 参考文献

1. COO格式规范: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
2. 稀疏矩阵算法: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems.
3. 裂变矩阵方法: Shaukat, N., et al. (2017). "Fission matrix approach for criticality calculations."
