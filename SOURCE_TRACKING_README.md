# OpenMC 源粒子标签和追踪系统使用说明

## 概述

本系统为OpenMC添加了源粒子标签和追踪功能，可以在活跃代中对抽样出来的源粒子打上标签、记录位置，并追踪其子代粒子。

## 功能特点

### 1. 源粒子标签
- **source_label**: 唯一标识符，格式为 `batch_id * 1000000 + particle_id`
- **source_position**: 记录源粒子的初始位置
- **source_batch**: 记录源粒子所属的batch号

### 2. 子代追踪
- 所有由源粒子产生的次级粒子（裂变、散射等）都会继承父粒子的源标签信息
- 可以追溯任意粒子的源头

### 3. 自动输出
- 在活跃代自动生成追踪文件 `source_tracking_batch_X.txt`
- 记录粒子的完整历史信息

## 代码修改内容

### 1. 数据结构修改

#### SourceSite 结构体（particle_data.h）
```cpp
struct SourceSite {
  // ...existing fields...
  
  // Custom labeling system for source particle tracking
  int64_t source_label {0};         //!< Custom label for source particle
  Position source_position;         //!< Initial source position
  int source_batch {-1};            //!< Batch number when source was created
};
```

#### ParticleData 类（particle_data.h）
```cpp
class ParticleData : public GeometryState {
private:
  // Source particle tracking variables
  int64_t source_label_ {0};         //!< Custom label for source particle tracking
  Position source_position_;         //!< Initial source position
  int source_batch_ {-1};            //!< Batch number when source was created

public:
  // Source particle tracking accessors
  int64_t& source_label() { return source_label_; }
  const int64_t& source_label() const { return source_label_; }
  Position& source_position() { return source_position_; }
  const Position& source_position() const { return source_position_; }
  int& source_batch() { return source_batch_; }
  const int& source_batch() const { return source_batch_; }
};
```

### 2. 初始化源标签（simulation.cpp）
```cpp
void initialize_history(Particle& p, int64_t index_source)
{
  // ...existing code...
  
  // Set source particle tracking labels (only for active generations)
  if (settings::run_mode == RunMode::EIGENVALUE && 
      simulation::current_batch > settings::n_inactive) {
    // Create unique source label: batch_id * max_particles + particle_id
    p.source_label() = static_cast<int64_t>(simulation::current_batch) * 1000000 + p.id();
    p.source_position() = p.r();  // Record initial position
    p.source_batch() = simulation::current_batch;
  }
}
```

### 3. 子代标签传递（particle.cpp）
- 在 `create_secondary()` 和 `split()` 函数中自动传递源标签
- 在 `from_source()` 函数中读取源标签信息

### 4. 追踪输出系统

#### 新增文件
- `include/openmc/source_tracking.h`: 头文件
- `src/source_tracking.cpp`: 实现文件

#### 功能
- 自动初始化追踪系统（仅活跃代）
- 在粒子碰撞时记录追踪信息
- 输出格式：`source_label batch source_x source_y source_z current_x current_y current_z energy weight collisions is_source`

## 使用方法

### 1. 编译
```bash
cd build
mingw32-make -j24
```

### 2. 运行
正常运行OpenMC，系统会在活跃代自动：
- 为源粒子分配标签
- 追踪所有相关粒子
- 输出追踪文件

### 3. 输出文件
- `source_tracking_batch_X.txt`: 每个活跃batch的追踪数据
- 包含完整的粒子历史和位置信息

## 输出文件格式

```
# Source Particle Tracking Data
# Columns: source_label source_batch source_x source_y source_z current_x current_y current_z energy weight collisions is_source
1000001 1 0.0 0.0 0.0 1.2 0.5 0.3 2.0e6 0.8 3 0
1000001 1 0.0 0.0 0.0 2.1 1.2 0.8 1.5e6 0.7 5 0
...
```

其中：
- `source_label`: 源粒子唯一标识
- `source_batch`: 源粒子所属batch
- `source_x/y/z`: 源粒子初始位置
- `current_x/y/z`: 当前粒子位置
- `energy`: 当前粒子能量
- `weight`: 当前粒子权重
- `collisions`: 碰撞次数
- `is_source`: 是否为原始源粒子（1=是，0=否）

## 注意事项

1. **性能影响**: 追踪功能仅在活跃代启用，对非活跃代无影响
2. **内存使用**: 系统会为每个粒子存储额外的标签信息
3. **文件输出**: 每个活跃batch会生成一个追踪文件
4. **MPI兼容**: 系统已考虑MPI并行运行的情况

## 扩展可能

1. **选择性追踪**: 可以添加条件来只追踪特定区域或特定类型的源粒子
2. **实时分析**: 可以在追踪过程中进行实时数据分析
3. **可视化接口**: 可以添加接口用于可视化粒子轨迹
4. **统计分析**: 可以添加源粒子贡献度分析等高级功能

这个系统为OpenMC提供了强大的源粒子追踪能力，可以用于详细的物理分析和调试。
