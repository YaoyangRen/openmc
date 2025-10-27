# OpenMC 格林函数功能使用手册

## 目录

1. [功能概述](#功能概述)
2. [理论背景](#理论背景)
3. [数据结构](#数据结构)
4. [API 参考](#api-参考)
5. [使用方法](#使用方法)
6. [数据格式](#数据格式)
7. [示例代码](#示例代码)
8. [性能优化](#性能优化)
9. [故障排除](#故障排除)

---

## 功能概述

OpenMC 格林函数功能实现了每个源粒子的独立追踪和统计，能够：

- ✅ 为每个源粒子生成独立的格林函数矩阵
- ✅ 自动累积所有源粒子的总体贡献
- ✅ 支持多批次（batch）计算和数据持久化
- ✅ 提供 HDF5 格式的结构化数据输出
- ✅ 线程安全的并行累积
- ✅ 自动检测模型边界或手动指定

### 主要应用场景

- **源重要性分析**: 评估不同源位置对系统响应的贡献
- **敏感性研究**: 分析源分布变化对结果的影响
- **响应函数计算**: 构建空间相关的响应函数
- **不确定性量化**: 评估源项不确定性传播
- **探测器优化**: 优化探测器位置和响应

---

## 理论背景

### 格林函数定义

格林函数 G(r, r') 表示位置 r' 处的单位点源对位置 r 处物理量的响应：

```
G(r, r') = φ(r | S(r'))
```

其中：
- `φ(r)` 是位置 r 处的物理量（如通量、反应率等）
- `S(r')` 是位置 r' 处的单位点源

### 线性叠加原理

任意源分布 S(r') 产生的响应可以通过格林函数叠加得到：

```
φ(r) = ∫ G(r, r') S(r') dr'
```

### 数值实现

在 OpenMC 中，通过追踪每个源粒子及其产生的所有次级粒子来构建格林函数：

```
G_i(r) = Σ [权重 × 路径长度] (源粒子i及其子代)
```

---

## 数据结构

### C++ 类定义

```cpp
class GreenFunctionMesh {
public:
  // 构造函数
  explicit GreenFunctionMesh(
    double resolution,                           // 网格分辨率 (cm)
    int max_batches,                             // 最大批次数
    bool auto_bounds = true,                     // 自动获取边界
    const std::array<double, 3>& manual_lower = {0.0, 0.0, 0.0},
    const std::array<double, 3>& manual_upper = {10.0, 10.0, 10.0}
  );

  // 累积贡献
  void accumulate(
    const Position& r,           // 空间位置
    double contribution,         // 贡献值
    int64_t source_particle_id   // 源粒子ID
  );

  // 批次管理
  void start_new_batch(int batch_id);
  
  // 数据输出
  void finalize_greenfunction_mesh(const int batch_id);
  
  // 数据查询
  const vector<double>& get_particle_data(int64_t source_particle_id) const;
  const std::array<int, 3>& shape() const;
  const std::array<double, 3>& origin() const;
  double pitch() const;

private:
  // 每个源粒子的格林函数矩阵
  std::unordered_map<int64_t, vector<double>> particle_green_functions_;
  
  // 当前批次的临时数据
  std::unordered_map<int64_t, vector<double>> current_batch_particle_data_;
  
  // 累积的总格林函数
  vector<double> cumulative_data_;
  
  // 网格参数
  std::array<int, 3> shape_;       // 网格维度
  std::array<double, 3> origin_;   // 网格原点
  double pitch_;                   // 网格间距
  double inv_pitch_;               // 间距倒数（优化）
  
  // 批次信息
  int current_batch_id_;
  int max_batches_;
  size_t spatial_size_;
  
  // 线程同步
  mutable std::mutex data_mutex_;
};
```

### 内存布局

格林函数数据以一维数组形式存储：

```
index = ix + shape[0] * (iy + shape[1] * iz)
```

其中 `(ix, iy, iz)` 是三维网格索引。

**空间复杂度**:
- 每个源粒子: `shape[0] × shape[1] × shape[2] × 8 bytes`
- N 个源粒子总计: `N × spatial_size × 8 bytes`

---

## API 参考

### 构造函数

```cpp
GreenFunctionMesh::GreenFunctionMesh(
  double resolution,
  int max_batches,
  bool auto_bounds = true,
  const std::array<double, 3>& manual_lower = {0.0, 0.0, 0.0},
  const std::array<double, 3>& manual_upper = {10.0, 10.0, 10.0}
)
```

**参数说明**:

| 参数 | 类型 | 说明 | 默认值 |
|------|------|------|--------|
| `resolution` | double | 网格单元边长 (cm) | 必填 |
| `max_batches` | int | 最大批次数量 | 必填 |
| `auto_bounds` | bool | 是否自动获取模型边界 | true |
| `manual_lower` | array<double,3> | 手动指定的下边界 | {0,0,0} |
| `manual_upper` | array<double,3> | 手动指定的上边界 | {10,10,10} |

**自动边界检测机制**:

当 `auto_bounds = true` 时，构造函数会：
1. 从 OpenMC 根宇宙（root universe）自动读取几何边界
2. 检查是否存在无限边界（如无限大平面）
3. 如果存在无限边界，回退到使用 `manual_lower` 和 `manual_upper`
4. 否则，在自动检测的边界上添加 `0.1 × pitch` 的边距，防止边缘粒子丢失

**边界来源**: 
- 几何信息从 XML 输入文件（`geometry.xml`）读取
- 通过 OpenMC 内部的 `model::universes` 和 `model::root_universe` 访问
- 对应于模型中所有 Cell 和 Surface 定义的包络边界

**示例**:

```cpp
// 自动获取边界，网格间距 1 cm（推荐）
auto gf_mesh = std::make_unique<GreenFunctionMesh>(1, 100);

// 手动指定边界
auto gf_mesh = std::make_unique<GreenFunctionMesh>(
  0.5, 100, false, 
  {-5.0, -5.0, -5.0},   // 下边界
  {15.0, 15.0, 15.0}    // 上边界
);
```

### accumulate 方法

```cpp
void accumulate(
  const Position& r,
  double contribution,
  int64_t source_particle_id
)
```

**功能**: 为特定源粒子在指定位置累积贡献值

**参数**:
- `r`: 粒子当前位置 (cm)
- `contribution`: 贡献值（通常为 权重 × 路径长度）
- `source_particle_id`: 源粒子唯一标识符

**线程安全**: 是（使用互斥锁和原子操作）

**边界处理**: 自动过滤网格范围外的贡献

### start_new_batch 方法

```cpp
void start_new_batch(int batch_id)
```

**功能**: 开始新的批次计算，保存上一批次数据

**调用时机**: 每个批次开始前

### finalize_greenfunction_mesh 方法

```cpp
void finalize_greenfunction_mesh(const int batch_id)
```

**功能**: 完成计算并将数据写入 HDF5 文件

**输出文件**: `green_function_data.h5`

**调用时机**: 模拟结束时

---

## 使用方法

### 1. 初始化格林函数网格

在 `simulation.cpp` 中初始化：

```cpp
// 在 initialize_batch() 中创建格林函数网格
if (!simulation::green_function_mesh) {
  simulation::green_function_mesh = std::make_unique<GreenFunctionMesh>(
    1.0,    // 1 cm 网格分辨率
    100,    // 最大 100 个批次
    true    // 自动获取模型边界
  );
}

// 开始新批次
if (settings::clutch_on && simulation::green_function_mesh) {
  simulation::green_function_mesh->start_new_batch(simulation::current_batch);
}
```

### 2. 粒子追踪中累积数据

在 `particle.cpp` 的 `event_advance()` 中：

```cpp
void Particle::event_advance() {
  // ...现有代码...
  
  // 移动粒子
  this->move_distance(distance);
  
  // 累积格林函数贡献
  if (simulation::green_function_mesh && material() != MATERIAL_VOID) {
    double contribution = wgt() * distance;  // 通量贡献
    simulation::green_function_mesh->accumulate(
      r(), contribution, source_particle_id()
    );
  }
  
  // ...现有代码...
}
```

### 3. 源粒子 ID 追踪

#### 设置源粒子 ID

在 `initialize_history()` 中：

```cpp
void initialize_history(Particle& p, int64_t index_source) {
  // ...现有代码...
  
  // 设置源粒子ID为当前粒子ID
  p.source_particle_id() = p.id();
  
  // ...现有代码...
}
```

#### 次级粒子 ID 继承

在 `create_secondary()` 中：

```cpp
bool Particle::create_secondary(double wgt, Direction u, double E, ParticleType type) {
  auto& bank = secondary_bank().emplace_back();
  
  // 复制物理状态
  bank.particle = type;
  bank.wgt = wgt;
  bank.r = r();
  bank.u = u;
  bank.E = E;
  
  // 继承源追踪信息
  bank.source_label = this->source_label();
  bank.source_position = this->source_position();
  bank.source_batch = this->source_batch();
  bank.source_particle_id = this->source_particle_id();  // 关键！
  
  return true;
}
```

### 4. 数据输出

在模拟结束时：

```cpp
// 在 finalize_batch() 或模拟结束时
if (settings::clutch_on && simulation::green_function_mesh) {
  if (simulation::current_batch == settings::n_batches) {
    simulation::green_function_mesh->finalize_greenfunction_mesh(
      simulation::current_batch
    );
  }
}
```

---

## 数据格式

### HDF5 文件结构

```
green_function_data.h5
├── 属性 (Attributes)
│   ├── filetype: "green_function_mesh_per_particle"
│   ├── version: "1.0"
│   ├── pitch: 1.0
│   └── n_source_particles: 10000
│
├── 数据集 (Datasets)
│   ├── shape [3]                      // 网格维度 [nx, ny, nz]
│   ├── origin [3]                     // 网格原点 [x0, y0, z0]
│   ├── cumulative_green_function [N]  // 累积格林函数
│   └── source_particle_ids [M]        // 源粒子ID列表
│
└── 群组 (Group)
    └── source_particles/
        ├── particle_1 [N]             // 源粒子1的格林函数
        ├── particle_10 [N]            // 源粒子10的格林函数
        └── ...                        // 其他源粒子
```

### 数据类型

| 数据集 | 类型 | 形状 | 说明 |
|--------|------|------|------|
| shape | int32 | (3,) | 网格在 X, Y, Z 方向的单元数 |
| origin | float64 | (3,) | 网格原点坐标 (cm) |
| cumulative_green_function | float64 | (N,) | 所有源粒子贡献的总和 |
| source_particle_ids | int64 | (M,) | 参与计算的源粒子ID列表 |
| particle_XXXX | float64 | (N,) | 单个源粒子的格林函数 |

其中 N = nx × ny × nz (总网格单元数), M = 源粒子数量

### 空间索引映射

**一维到三维**:
```python
gf_3d = gf_1d.reshape(shape)
value = gf_3d[ix, iy, iz]
```

**三维到一维**:
```python
index = ix + shape[0] * (iy + shape[1] * iz)
value = gf_1d[index]
```

**物理坐标计算**:
```python
x = origin[0] + ix * pitch
y = origin[1] + iy * pitch
z = origin[2] + iz * pitch
```

---

## 示例代码

### Python 数据读取示例

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# 读取HDF5文件
with h5py.File('green_function_data.h5', 'r') as f:
    # 读取网格信息
    shape = f['shape'][:]
    origin = f['origin'][:]
    pitch = f.attrs['pitch']
    
    # 读取累积格林函数
    cumulative_gf = f['cumulative_green_function'][:]
    cumulative_3d = cumulative_gf.reshape(shape)
    
    # 读取源粒子ID列表
    particle_ids = f['source_particle_ids'][:]
    
    # 读取特定源粒子的数据
    pid = particle_ids[0]
    particle_gf = f[f'source_particles/particle_{pid}'][:]
    particle_3d = particle_gf.reshape(shape)
    
    # 验证线性叠加
    total_from_particles = 0.0
    for pid in particle_ids:
        data = f[f'source_particles/particle_{pid}'][:]
        total_from_particles += np.sum(data)
    
    print(f"累积总和: {np.sum(cumulative_gf):.6e}")
    print(f"单粒子总和: {total_from_particles:.6e}")
    print(f"相对误差: {abs(total_from_particles - np.sum(cumulative_gf)) / np.sum(cumulative_gf) * 100:.4f}%")
```

### 可视化示例

```python
# 绘制中心截面
center_z = shape[2] // 2
center_slice = cumulative_3d[:, :, center_z]

# 创建坐标网格
x = np.linspace(origin[0], origin[0] + shape[0] * pitch, shape[0])
y = np.linspace(origin[1], origin[1] + shape[1] * pitch, shape[1])
X, Y = np.meshgrid(x, y)

# 绘图
plt.figure(figsize=(10, 8))
plt.contourf(X, Y, center_slice.T, levels=50, cmap='viridis')
plt.colorbar(label='格林函数值')
plt.xlabel('X (cm)')
plt.ylabel('Y (cm)')
plt.title(f'累积格林函数中心截面 (z = {origin[2] + center_z * pitch:.2f} cm)')
plt.grid(True, alpha=0.3)
plt.savefig('green_function_slice.png', dpi=300)
plt.show()
```

### 源重要性分析示例

```python
# 分析每个源粒子的贡献
particle_contributions = {}

with h5py.File('green_function_data.h5', 'r') as f:
    particle_ids = f['source_particle_ids'][:]
    
    for pid in particle_ids:
        data = f[f'source_particles/particle_{pid}'][:]
        contribution = np.sum(data)
        particle_contributions[pid] = contribution

# 排序并显示前10个最重要的源粒子
sorted_particles = sorted(particle_contributions.items(), 
                         key=lambda x: x[1], reverse=True)

print("前10个最重要的源粒子:")
total_contribution = sum(particle_contributions.values())
for i, (pid, contrib) in enumerate(sorted_particles[:10]):
    percentage = contrib / total_contribution * 100
    print(f"{i+1}. 源粒子 {pid}: {contrib:.6e} ({percentage:.2f}%)")
```

---

## 性能优化

### 1. 内存管理

**稀疏存储**:
- 使用 `std::unordered_map` 仅存储有贡献的源粒子
- 自动过滤零贡献粒子

**批次处理**:
- 分批次累积数据，减少内存峰值
- 每批次结束时清空临时缓存

**内存估算**:
```
总内存 ≈ N_particles × (nx × ny × nz) × 8 bytes

例如: 10,000 粒子 × (100 × 100 × 100) × 8 bytes ≈ 76 GB
```

### 2. 并行性能

**线程安全机制**:
- 使用 `std::mutex` 保护 `unordered_map` 访问
- 使用 `#pragma omp atomic` 保护数值累积
- 细粒度锁，最小化锁持有时间

**并行策略**:
```cpp
{
  // 短暂持有锁获取数据指针
  std::lock_guard<std::mutex> lock(data_mutex_);
  auto it = current_batch_particle_data_.find(source_particle_id);
  particle_data_ptr = &(it->second);
}

// 在锁外进行原子累积
#pragma omp atomic
(*particle_data_ptr)[index] += contribution;
```

### 3. 计算优化

**预计算优化**:
- 存储 `inv_pitch_` = 1.0 / `pitch_`，避免重复除法
- 使用整数索引计算，避免浮点运算

**边界检查优化**:
```cpp
// 快速边界检查
if (ix >= 0 && ix < shape_[0] && 
    iy >= 0 && iy < shape_[1] && 
    iz >= 0 && iz < shape_[2]) {
  // 在范围内才进行索引计算
}
```

### 4. I/O 优化

**HDF5 写入**:
- 批量写入数据集，减少I/O调用
- 使用压缩（可选）减小文件大小

**文件大小估算**:
```
文件大小 ≈ (N_particles + 1) × spatial_size × 8 bytes + 元数据

例如: (10,000 + 1) × (100³) × 8 ≈ 76 GB
```

---

## 故障排除

### 常见问题

#### 1. "Particle data size mismatch" 警告

**原因**: 多线程并发访问导致数据结构不一致

**解决方案**: 
- 已通过 `std::mutex` 保护解决
- 确保使用最新版本的代码

#### 2. "Index out of bounds" 错误

**原因**: 粒子位置超出网格边界

**检查**:
```cpp
// 查看初始化输出
GreenFunctionMesh initialized:
  Bounds: [x_min, y_min, z_min] to [x_max, y_max, z_max]
```

**解决方案**:
- 增加网格边界扩展 (`margin`)
- 检查模型边界是否正确

#### 3. 内存不足

**症状**: 程序崩溃或运行缓慢

**解决方案**:
- 增大网格间距 (`pitch`)
- 减少源粒子数量
- 使用分批次处理
- 考虑使用稀疏存储格式

#### 4. 文件无法生成

**检查**:
```cpp
// 确保调用了finalize
simulation::green_function_mesh->finalize_greenfunction_mesh(batch_id);
```

**权限**: 确保程序有写入权限

### 调试技巧

#### 启用详细输出

```cpp
// 在构造函数中已有调试信息
std::cout << "GreenFunctionMesh initialized:" << std::endl;
std::cout << "  Bounds: ..." << std::endl;
std::cout << "  Shape: ..." << std::endl;
```

#### 验证数据一致性

```python
# Python验证脚本
import h5py
import numpy as np

with h5py.File('green_function_data.h5', 'r') as f:
    cumulative = f['cumulative_green_function'][:]
    particle_ids = f['source_particle_ids'][:]
    
    total = 0.0
    for pid in particle_ids:
        data = f[f'source_particles/particle_{pid}'][:]
        total += np.sum(data)
    
    diff = abs(total - np.sum(cumulative))
    rel_error = diff / np.sum(cumulative) * 100
    
    print(f"线性叠加验证:")
    print(f"  累积总和: {np.sum(cumulative):.6e}")
    print(f"  单粒子总和: {total:.6e}")
    print(f"  差值: {diff:.6e}")
    print(f"  相对误差: {rel_error:.4f}%")
    
    if rel_error < 0.01:
        print("  ✓ 验证通过")
    else:
        print("  ✗ 验证失败")
```

---

## 附录

### A. 完整的 C++ 接口

```cpp
namespace openmc {

class GreenFunctionMesh {
public:
  explicit GreenFunctionMesh(
    double resolution,
    int max_batches,
    bool auto_bounds = true,
    const std::array<double, 3>& manual_lower = {0.0, 0.0, 0.0},
    const std::array<double, 3>& manual_upper = {10.0, 10.0, 10.0}
  );

  void accumulate(
    const Position& r,
    double contribution,
    int64_t source_particle_id
  );

  void start_new_batch(int batch_id);
  void finalize_greenfunction_mesh(const int batch_id);
  
  const vector<double>& get_particle_data(int64_t source_particle_id) const;
  const std::array<int, 3>& shape() const { return shape_; }
  const std::array<double, 3>& origin() const { return origin_; }
  double pitch() const { return pitch_; }
};

// 全局实例
namespace simulation {
  extern std::unique_ptr<GreenFunctionMesh> green_function_mesh;
}

} // namespace openmc
```

### B. Python 分析工具库

建议创建一个 Python 模块用于格林函数数据分析：

```python
# greenfunction_analysis.py

import h5py
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional

class GreenFunctionAnalyzer:
    """格林函数数据分析器"""
    
    def __init__(self, filename: str = 'green_function_data.h5'):
        self.filename = filename
        self._load_metadata()
    
    def _load_metadata(self):
        """加载元数据"""
        with h5py.File(self.filename, 'r') as f:
            self.shape = f['shape'][:]
            self.origin = f['origin'][:]
            self.pitch = f.attrs['pitch']
            self.n_particles = f.attrs['n_source_particles']
            self.particle_ids = f['source_particle_ids'][:]
    
    def get_cumulative(self) -> np.ndarray:
        """获取累积格林函数"""
        with h5py.File(self.filename, 'r') as f:
            data = f['cumulative_green_function'][:]
        return data.reshape(self.shape)
    
    def get_particle_data(self, particle_id: int) -> np.ndarray:
        """获取特定源粒子的格林函数"""
        with h5py.File(self.filename, 'r') as f:
            data = f[f'source_particles/particle_{particle_id}'][:]
        return data.reshape(self.shape)
    
    def verify_superposition(self, sample_size: int = 100) -> dict:
        """验证线性叠加性质"""
        with h5py.File(self.filename, 'r') as f:
            cumulative_sum = np.sum(f['cumulative_green_function'][:])
            
            # 采样验证
            sample_ids = np.random.choice(
                self.particle_ids, 
                min(sample_size, len(self.particle_ids)), 
                replace=False
            )
            
            particle_sum = 0.0
            for pid in sample_ids:
                particle_sum += np.sum(f[f'source_particles/particle_{pid}'][:])
            
            # 估算总和
            estimated_total = particle_sum * len(self.particle_ids) / len(sample_ids)
        
        return {
            'cumulative_sum': cumulative_sum,
            'estimated_sum': estimated_total,
            'difference': abs(cumulative_sum - estimated_total),
            'relative_error': abs(cumulative_sum - estimated_total) / cumulative_sum * 100
        }
    
    def plot_slice(self, axis: str = 'z', index: Optional[int] = None, 
                   figsize: Tuple[int, int] = (10, 8)):
        """绘制截面图"""
        data_3d = self.get_cumulative()
        
        if index is None:
            index = self.shape[{'x': 0, 'y': 1, 'z': 2}[axis]] // 2
        
        if axis == 'z':
            slice_data = data_3d[:, :, index]
            extent = [self.origin[0], self.origin[0] + self.shape[0] * self.pitch,
                     self.origin[1], self.origin[1] + self.shape[1] * self.pitch]
            xlabel, ylabel = 'X (cm)', 'Y (cm)'
        elif axis == 'y':
            slice_data = data_3d[:, index, :]
            extent = [self.origin[0], self.origin[0] + self.shape[0] * self.pitch,
                     self.origin[2], self.origin[2] + self.shape[2] * self.pitch]
            xlabel, ylabel = 'X (cm)', 'Z (cm)'
        else:  # x
            slice_data = data_3d[index, :, :]
            extent = [self.origin[1], self.origin[1] + self.shape[1] * self.pitch,
                     self.origin[2], self.origin[2] + self.shape[2] * self.pitch]
            xlabel, ylabel = 'Y (cm)', 'Z (cm)'
        
        plt.figure(figsize=figsize)
        plt.imshow(slice_data.T, extent=extent, origin='lower', 
                  cmap='viridis', aspect='auto')
        plt.colorbar(label='格林函数值')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(f'累积格林函数 {axis.upper()} 截面 (index={index})')
        plt.grid(True, alpha=0.3)
        return plt.gcf()
```

### C. 参考文献

1. Bell, G. I., & Glasstone, S. (1970). *Nuclear Reactor Theory*. Van Nostrand Reinhold.

2. Duderstadt, J. J., & Hamilton, L. J. (1976). *Nuclear Reactor Analysis*. John Wiley & Sons.

3. Romano, P. K., et al. (2015). "OpenMC: A state-of-the-art Monte Carlo code for research and development." *Annals of Nuclear Energy*, 82, 90-97.

4. Williams, M. L. (1986). "Generalized Contributon Response Theory." *Nuclear Science and Engineering*, 93(3), 288-297.

---

## 版本历史

| 版本 | 日期 | 更新内容 |
|------|------|----------|
| 1.0 | 2025-10-09 | 初始版本，实现基本格林函数功能 |
| 1.1 | 2025-10-09 | 添加自动边界检测和线程安全优化 |

---

**最后更新**: 2025年10月9日
