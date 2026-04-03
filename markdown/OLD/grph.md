# 传递函数与裂变矩阵算法流程图

## 1. 总体流程概览

```mermaid
flowchart TD
    Start([开始模拟]) --> Init[初始化阶段]
    Init --> Batch[批次循环]
    Batch --> Particle[粒子模拟循环]
    Particle --> Finalize[最终化处理]
    Finalize --> End([结束])
    
    Init -.->|详见图2| InitDetail
    Particle -.->|详见图3| ParticleDetail
    Particle -.->|格林函数详见图4| GFDetail
    Particle -.->|裂变矩阵详见图5| FMDetail
    Finalize -.->|详见图6| FinalDetail
    
    style Start fill:#90EE90
    style End fill:#FFB6C1
    style Init fill:#FFE4B5
    style Particle fill:#87CEEB
    style Finalize fill:#DDA0DD
```

---

## 2. 初始化阶段详细流程

```mermaid
flowchart TD
    Start([初始化开始]) --> LoadModel[加载几何模型]
    LoadModel --> LoadMaterial[加载材料数据]
    LoadMaterial --> InitSettings[读取模拟设置]
    
    InitSettings --> CheckGF{启用格林函数?}
    CheckGF -->|是| CreateGF[创建GreenFunctionMesh]
    
    CreateGF --> GetBounds[从root_universe获取边界]
    GetBounds --> CalcShape[计算网格shape<br/>shape = ceil边界-原点/pitch]
    CalcShape --> AllocGF[初始化数据结构<br/>- particle_green_functions_<br/>- current_batch_data_<br/>- cumulative_data_]
    
    CheckGF -->|否| CheckFM
    AllocGF --> CheckFM
    
    CheckFM{启用裂变矩阵?} -->|是| CreateFM[创建FissionMatrix]
    
    CreateFM --> CalcFMShape[计算网格维度<br/>n_cells = nx × ny × nz]
    CalcFMShape --> AllocFM[初始化稀疏矩阵<br/>- fission_matrix_sparse_<br/>- current_batch_sparse_<br/>- source_birth_cells_<br/>- source_counts_]
    
    CheckFM -->|否| InitComplete
    AllocFM --> InitComplete[初始化完成]
    
    InitComplete --> Return([返回主循环])
    
    style Start fill:#90EE90
    style Return fill:#87CEEB
    style CreateGF fill:#87CEEB
    style CreateFM fill:#DDA0DD
```

---

## 3. 单粒子输运主循环

```mermaid
flowchart TD
    Start([开始粒子历史]) --> InitParticle[初始化粒子<br/>from_source]
    
    InitParticle --> SetSourceID[设置source_particle_id]
    SetSourceID --> RecordBirth{裂变矩阵启用?}
    
    RecordBirth -->|是| RecordSource[记录源粒子出生<br/>source_birth_cells_<br/>source_id → cell_index]
    RecordBirth -->|否| Transport
    RecordSource --> Transport
    
    Transport[传输循环开始] --> CalcXS[计算宏观截面<br/>event_calculate_xs]
    CalcXS --> SampleDistance[采样自由程<br/>-ln ξ/Σt]
    
    SampleDistance --> CalcBoundary[计算到边界距离]
    CalcBoundary --> MinDistance[取最小距离<br/>distance = min]
    
    MinDistance --> MoveParticle[移动粒子<br/>r = r + u×distance]
    MoveParticle --> AccumGF{格林函数启用?}
    
    AccumGF -->|是| CallGFAccum[调用格林函数累积<br/>详见图4]
    AccumGF -->|否| CheckEvent
    CallGFAccum --> CheckEvent
    
    CheckEvent{事件类型?}
    CheckEvent -->|边界| HandleBoundary[边界处理<br/>详见图3.1]
    CheckEvent -->|碰撞| HandleCollision[碰撞处理<br/>详见图3.2]
    
    HandleBoundary --> CheckAlive
    HandleCollision --> CheckAlive
    
    CheckAlive{粒子存活?}
    CheckAlive -->|是| Transport
    CheckAlive -->|否| CheckSecondary{有次级粒子?}
    
    CheckSecondary -->|是| ReviveSecondary[复活次级粒子<br/>继承source_particle_id]
    CheckSecondary -->|否| End([粒子历史结束])
    
    ReviveSecondary --> Transport
    
    style Start fill:#90EE90
    style End fill:#FFB6C1
    style Transport fill:#87CEEB
    style AccumGF fill:#FFA07A
```

### 3.1 边界处理子流程

```mermaid
flowchart TD
    Start([到达边界]) --> CrossSurface[event_cross_surface]
    CrossSurface --> CheckBC{边界条件?}
    
    CheckBC -->|真空| VacuumBC[粒子死亡<br/>wgt = 0]
    CheckBC -->|反射| ReflectBC[反射粒子<br/>u = u - 2n·u n]
    CheckBC -->|周期| PeriodicBC[周期边界<br/>r = r + Δr]
    CheckBC -->|透射| TransBC[继续传输]
    
    VacuumBC --> Return([返回])
    ReflectBC --> Return
    PeriodicBC --> Return
    TransBC --> Return
    
    style Start fill:#FFE4B5
    style Return fill:#87CEEB
```

### 3.2 碰撞处理子流程

```mermaid
flowchart TD
    Start([发生碰撞]) --> SampleReaction[采样反应类型<br/>根据Σs, Σf, Σa]
    
    SampleReaction --> ReactionType{反应类型?}
    
    ReactionType -->|散射| Scatter[散射处理]
    ReactionType -->|吸收| Absorb[吸收<br/>wgt = 0]
    ReactionType -->|裂变| Fission[裂变处理]
    
    Scatter --> SampleEnergy[采样出射能量<br/>根据散射律]
    SampleEnergy --> SampleAngle[采样散射角度]
    SampleAngle --> UpdateDirection[更新方向u]
    UpdateDirection --> Return
    
    Absorb --> Return([返回])
    
    Fission --> CalcNu[计算ν裂变中子数]
    CalcNu --> RecordFM{裂变矩阵启用?}
    
    RecordFM -->|是| CallFMRecord[调用裂变矩阵记录<br/>详见图5]
    RecordFM -->|否| CreateSecondary
    CallFMRecord --> CreateSecondary
    
    CreateSecondary --> BankLoop{遍历ν个次级粒子}
    BankLoop -->|下一个| SampleFission[采样裂变中子<br/>能量、方向]
    SampleFission --> InheritSource[继承源信息<br/>source_particle_id<br/>source_position]
    InheritSource --> BankSecondary[存入secondary_bank]
    BankSecondary --> BankLoop
    
    BankLoop -->|完成| Return
    
    style Start fill:#FFE4B5
    style Return fill:#87CEEB
    style Fission fill:#DDA0DD
```

---

## 4. 格林函数累积详细流程

```mermaid
flowchart TD
    Start([格林函数累积]) --> CalcContrib[计算贡献<br/>contribution = w×distance×νΣf/Σt]
    
    CalcContrib --> CheckPositive{contribution > 0?}
    CheckPositive -->|否| Return([返回])
    
    CheckPositive -->|是| PosToIndex[位置转索引<br/>ix = floor x-x0/pitch<br/>iy = floor y-y0/pitch<br/>iz = floor z-z0/pitch]
    
    PosToIndex --> CalcLinear[计算线性索引<br/>index = ix + nx×iy + nx×ny×iz]
    
    CalcLinear --> BoundsCheck{在网格内?}
    BoundsCheck -->|否| DropContrib[丢弃贡献<br/>dropped_contributions_++]
    DropContrib --> Return
    
    BoundsCheck -->|是| LockMutex[加锁data_mutex_]
    
    LockMutex --> FindData{粒子数据存在?<br/>source_particle_id}
    
    FindData -->|否| CreateData[创建新数据<br/>current_batch_data_<br/>source_id → vector spatial_size_]
    FindData -->|是| GetData[获取现有数据]
    
    CreateData --> Unlock[解锁]
    GetData --> Unlock
    
    Unlock --> AtomicAdd1[原子累积<br/>particle_data index<br/>+= contribution]
    
    AtomicAdd1 --> AtomicAdd2[原子累积<br/>cumulative_data index<br/>+= contribution]
    
    AtomicAdd2 --> IncCounter[total_contributions_++]
    IncCounter --> Return
    
    style Start fill:#87CEEB
    style Return fill:#87CEEB
    style AtomicAdd1 fill:#FFA07A
    style AtomicAdd2 fill:#FFA07A
    style LockMutex fill:#F0E68C
```

---

## 5. 裂变矩阵记录详细流程

```mermaid
flowchart TD
    Start([裂变矩阵记录]) --> GetFissionPos[获取裂变位置<br/>r_fission = particle.r]
    
    GetFissionPos --> FissionToIndex[裂变位置转索引<br/>j_fission = position_to_index r]
    
    FissionToIndex --> CheckFissionValid{j_fission有效?}
    CheckFissionValid -->|否| Return([返回])
    
    CheckFissionValid -->|是| LockMutex[加锁data_mutex_]
    
    LockMutex --> FindSource{查找源位置<br/>source_birth_cells_<br/>source_particle_id}
    
    FindSource -->|未找到| Unlock1[解锁]
    Unlock1 --> Return
    
    FindSource -->|找到| GetSourceCell[i_source = cell_index]
    
    GetSourceCell --> CalcKey[计算稀疏矩阵键<br/>key = i_source × n_cells + j_fission]
    
    CalcKey --> UpdateSparse[更新稀疏矩阵<br/>current_batch_sparse_<br/>key += nu_fission]
    
    UpdateSparse --> IncFissions[total_fissions_++]
    
    IncFissions --> Unlock2[解锁]
    Unlock2 --> Return
    
    style Start fill:#DDA0DD
    style Return fill:#DDA0DD
    style UpdateSparse fill:#FFA07A
    style LockMutex fill:#F0E68C
```

---

## 6. 最终化处理详细流程

```mermaid
flowchart TD
    Start([最终化开始]) --> SaveLastBatch[保存最后batch<br/>start_new_batch-1]
    
    SaveLastBatch --> ProcessGF{处理格林函数?}
    
    ProcessGF -->|是| OpenGFFile[创建HDF5文件<br/>transfer_function_data.h5]
    
    OpenGFFile --> WriteGFMeta[写入元数据<br/>- filetype: transfer_function<br/>- version: 2.0<br/>- pitch, origin, shape]
    
    WriteGFMeta --> WriteCumulative[写入累积传递函数<br/>cumulative_transfer_function]
    
    WriteCumulative --> CreateParticlesGroup[创建source_particles组]
    
    CreateParticlesGroup --> ParticleLoop{遍历源粒子}
    ParticleLoop -->|下一个| WriteParticleData[写入粒子数据<br/>particle_source_id]
    WriteParticleData --> ParticleLoop
    
    ParticleLoop -->|完成| WriteIDs[写入粒子ID列表<br/>source_particle_ids]
    WriteIDs --> CloseGF[关闭HDF5文件]
    
    ProcessGF -->|否| ProcessFM
    CloseGF --> ProcessFM
    
    ProcessFM{处理裂变矩阵?}
    
    ProcessFM -->|是| NormalizeFM[归一化裂变矩阵<br/>F_ij /= source_counts_i]
    
    NormalizeFM --> ConvertCOO[转换为COO格式]
    
    ConvertCOO --> COOLoop{遍历稀疏元素}
    COOLoop --> ExtractIndices[提取行列索引<br/>i = key / n_cells<br/>j = key % n_cells]
    ExtractIndices --> StoreEntry[存储<br/>row_indices i<br/>col_indices j<br/>data_raw value<br/>data_normalized value/count]
    StoreEntry --> COOLoop
    
    COOLoop -->|完成| OpenFMFile[创建HDF5文件<br/>fission_matrix.h5]
    
    OpenFMFile --> WriteFMMeta[写入元数据<br/>- storage_format: COO<br/>- nnz, sparsity<br/>- n_cells, pitch]
    
    WriteFMMeta --> WriteCOO[写入COO数据<br/>- row_indices<br/>- col_indices<br/>- data_raw<br/>- data_normalized]
    
    WriteCOO --> WriteSourceCounts[写入source_counts]
    
    WriteSourceCounts --> CalcStats[计算统计信息<br/>- 稀疏度<br/>- 内存节省<br/>- 最大元素]
    
    CalcStats --> CloseFM[关闭HDF5文件]
    
    ProcessFM -->|否| PrintSummary
    CloseFM --> PrintSummary
    
    PrintSummary[输出总结信息] --> End([最终化完成])
    
    style Start fill:#DDA0DD
    style End fill:#FFB6C1
    style ProcessGF fill:#87CEEB
    style ProcessFM fill:#DDA0DD
    style ConvertCOO fill:#98FB98
```

---

## 7. 批次管理流程

```mermaid
flowchart TD
    Start([开始新批次]) --> CheckPrevious{上一批次ID ≥ 0?}
    
    CheckPrevious -->|否| ClearCurrent[清空current数据]
    
    CheckPrevious -->|是| MergeGF{格林函数启用?}
    
    MergeGF -->|是| LoopParticles{遍历current_batch_data_}
    
    LoopParticles --> CheckTotal[计算粒子总贡献<br/>sum particle_data]
    CheckTotal --> HasData{总贡献 > 0?}
    
    HasData -->|否| LoopParticles
    
    HasData -->|是| CreatePerm{永久数据存在?}
    CreatePerm -->|否| InitPerm[创建particle_green_functions_<br/>source_id → zeros]
    CreatePerm -->|是| GetPerm[获取现有数据]
    
    InitPerm --> Accumulate[累积到永久存储<br/>permanent += current]
    GetPerm --> Accumulate
    
    Accumulate --> LoopParticles
    LoopParticles -->|完成| MergeFM
    
    MergeGF -->|否| MergeFM
    
    MergeFM{裂变矩阵启用?}
    
    MergeFM -->|是| LoopSparse{遍历current_batch_sparse_}
    
    LoopSparse --> MergeEntry[fission_matrix_sparse_<br/>key += value]
    MergeEntry --> LoopSparse
    
    LoopSparse -->|完成| MergeSource[合并source_counts_<br/>counts += current_counts]
    
    MergeSource --> IncRealization[n_realizations_++]
    
    MergeFM -->|否| ClearCurrent
    IncRealization --> ClearCurrent
    
    ClearCurrent --> ClearData[清空current_batch_data_<br/>清空current_batch_sparse_<br/>清空source_birth_cells_]
    
    ClearData --> SetBatchID[current_batch_id_ = new_id]
    SetBatchID --> End([批次切换完成])
    
    style Start fill:#FFE4B5
    style End fill:#87CEEB
    style MergeGF fill:#87CEEB
    style MergeFM fill:#DDA0DD
```

---

## 8. 数据结构关系图

```mermaid
flowchart LR
    subgraph GFData["格林函数数据结构"]
        GF1["particle_green_functions_
        map: source_id → vector(double)"]
        GF2["current_batch_particle_data_
        map: source_id → vector(double)"]
        GF3["cumulative_data_
        vector(double): spatial_size_"]
    end

    subgraph FMData["裂变矩阵数据结构 - 稀疏存储"]
        FM1["fission_matrix_sparse_
        map: key → value
        key = i × n + j"]
        FM2["current_batch_sparse_
        map: key → value"]
        FM3["source_birth_cells_
        map: source_id → cell_index"]
        FM4["source_counts_
        vector(double): n_cells"]
    end
    
    subgraph ParticleData["粒子数据"]
        P1["Particle
        - source_particle_id_
        - r_ (position)
        - wgt_ (weight)
        - E_ (energy)"]
    end
    
    P1 -->|传输过程累积| GF2
    P1 -->|出生时记录| FM3
    P1 -->|裂变时累积| FM2
    
    GF2 -->|batch结束合并| GF1
    GF1 -->|最终化累加| GF3
    
    FM2 -->|batch结束合并| FM1
    FM3 -->|查询源位置| FM1
    FM1 -->|归一化计算| FM4
    
    style GFData fill:#E6F3FF
    style FMData fill:#FFE6F0
    style ParticleData fill:#E6FFE6
```

---

## 9. 时序图：单粒子历史的格林函数与裂变矩阵交互

```mermaid
sequenceDiagram
    participant Main as 主循环
    participant Particle as 粒子
    participant GF as GreenFunctionMesh
    participant FM as FissionMatrix
    
    Main->>Particle: initialize_history(id)
    Particle->>Particle: source_particle_id = id
    
    alt 裂变矩阵启用
        Particle->>FM: record_source_birth(r, source_id)
        FM->>FM: source_birth_cells_[source_id] = cell_index
        FM->>FM: source_counts[cell_index] += 1
    end
    
    loop 传输循环
        Particle->>Particle: event_advance()
        Particle->>Particle: move(distance)
        
        alt 格林函数启用
            Particle->>GF: accumulate(r, contribution, source_id)
            GF->>GF: 加锁mutex
            GF->>GF: 获取/创建 particle_data[source_id]
            GF->>GF: 解锁mutex
            GF->>GF: atomic: particle_data[index] += contribution
            GF->>GF: atomic: cumulative_data[index] += contribution
        end
        
        alt 发生碰撞
            Particle->>Particle: event_collide()
            
            alt 裂变反应
                Particle->>Particle: 计算 nu_fission
                
                alt 裂变矩阵启用
                    Particle->>FM: record_fission_event(r, nu_fission, source_id)
                    FM->>FM: 加锁mutex
                    FM->>FM: source_cell = source_birth_cells_[source_id]
                    FM->>FM: fission_cell = position_to_index(r)
                    FM->>FM: key = source_cell × n + fission_cell
                    FM->>FM: current_batch_sparse_[key] += nu_fission
                    FM->>FM: 解锁mutex
                end
                
                Particle->>Particle: create_secondary()
                Particle->>Particle: 次级粒子继承 source_particle_id
            end
        end
        
        alt 粒子死亡且有次级粒子
            Particle->>Particle: revive_from_secondary()
            Particle->>Particle: 保持相同 source_particle_id
        end
    end
    
    Particle-->>Main: 粒子历史结束
```

---

## 10. 数据处理流程总览

```mermaid
flowchart TB
    subgraph Input["输入阶段"]
        I1[源粒子参数
        位置、能量、方向]
        I2[几何模型
        材料、边界]
        I3[网格参数
        分辨率、边界]
    end
    
    subgraph Processing["粒子模拟阶段"]
        P1[粒子传输
        位置变化]
        P2[碰撞处理
        反应采样]
        P3[次级粒子
        产生与追踪]
        
        P1 --> P2
        P2 --> P3
        P3 --> P1
    end
    
    subgraph Accumulation["累积阶段"]
        A1[格林函数累积
        T(source_id, r)]
        A2[裂变矩阵累积
        F(i, j)]
        A3[批次管理
        batch统计]
        
        A1 --> A3
        A2 --> A3
    end
    
    subgraph Output["输出阶段"]
        O1[归一化处理
        按源粒子数/批次数]
        O2[格式转换
        COO稀疏格式]
        O3[HDF5文件输出
        可视化数据]
    end
    
    I1 --> Processing
    I2 --> Processing
    I3 --> Accumulation
    
    Processing -->|位置+权重+νΣf| A1
    Processing -->|源位置+裂变位置| A2
    
    Accumulation --> O1
    O1 --> O2
    O2 --> O3
    
    style Input fill:#E6F3FF
    style Processing fill:#FFF4E6
    style Accumulation fill:#E6FFE6
    style Output fill:#FFE6F0
```

---

## 附录：核心算法伪代码

### A1. 格林函数累积

```python
function accumulate_green_function(particle):
    """
    在每次粒子推进时调用
    计算传递函数 T(P0 -> r) = ∫∫ ν̄Σf Φ dΩdE
    """
    # 计算贡献：w × distance × (ν̄Σf / Σt)
    contribution = particle.weight × distance × (nu_sigma_f / sigma_t)
    
    if contribution <= 0:
        return
    
    # 位置转换为网格索引
    cell_index = position_to_index(particle.position)
    
    if cell_index < 0:  # 超出边界
        dropped_contributions++
        return
    
    # 线程安全地获取或创建粒子数据
    with mutex_lock:
        if particle.source_id not in current_batch_data:
            current_batch_data[particle.source_id] = zeros(spatial_size)
    
    # 原子操作累积
    atomic_add(current_batch_data[particle.source_id][cell_index], contribution)
    atomic_add(cumulative_data[cell_index], contribution)
    total_contributions++
```

### A2. 裂变矩阵记录

```python
function record_source_birth(particle):
    """粒子初始化时记录源位置"""
    source_cell = position_to_index(particle.birth_position)
    
    with mutex_lock:
        source_birth_cells[particle.source_id] = source_cell
        source_counts[source_cell] += 1
        total_sources++

function record_fission_event(particle):
    """裂变事件时记录"""
    # 查找源位置
    with mutex_lock:
        if particle.source_id not in source_birth_cells:
            return
        source_cell = source_birth_cells[particle.source_id]
    
    # 计算裂变位置
    fission_cell = position_to_index(particle.position)
    if fission_cell < 0:
        return
    
    # 计算裂变产生中子数
    nu_fission = calculate_nu_fission(particle)
    
    # 累积到稀疏裂变矩阵
    key = source_cell × n_cells + fission_cell
    
    with mutex_lock:
        current_batch_sparse[key] += nu_fission
        total_fissions++
```

### A3. 批次管理

```python
function start_new_batch(batch_id):
    """批次切换时合并数据"""
    if current_batch_id >= 0:
        # 格林函数：保存上一个batch的数据
        for source_id, data in current_batch_data:
            if sum(data) > 0:
                if source_id not in particle_green_functions:
                    particle_green_functions[source_id] = zeros(spatial_size)
                # 累积到永久存储
                particle_green_functions[source_id] += data
        
        # 裂变矩阵：合并稀疏数据
        for key, value in current_batch_sparse:
            fission_matrix_sparse[key] += value
        
        for i in range(n_cells):
            source_counts[i] += current_batch_source_counts[i]
        
        n_realizations++
    
    # 清空当前batch数据
    current_batch_data.clear()
    current_batch_sparse.clear()
    source_birth_cells.clear()
    current_batch_source_counts.fill(0)
    
    current_batch_id = batch_id
```

### A4. COO格式转换

```python
function convert_to_COO(fission_matrix_sparse):
    """
    将稀疏矩阵转换为COO (Coordinate) 格式
    COO格式: (row_indices[], col_indices[], data[])
    """
    nnz = len(fission_matrix_sparse)
    
    row_indices = zeros(nnz, int)
    col_indices = zeros(nnz, int)
    data_raw = zeros(nnz, double)
    data_normalized = zeros(nnz, double)
    
    index = 0
    for key, value in fission_matrix_sparse:
        # 从线性键恢复(i,j)索引
        i = key // n_cells  # 行索引（源单元）
        j = key % n_cells   # 列索引（裂变单元）
        
        row_indices[index] = i
        col_indices[index] = j
        data_raw[index] = value
        
        # 归一化：F_norm[i][j] = F_raw[i][j] / source_counts[i]
        if source_counts[i] > 0:
            data_normalized[index] = value / source_counts[i]
        else:
            data_normalized[index] = 0
        
        index++
    
    return (row_indices, col_indices, data_raw, data_normalized)
```

---

## 流程图使用说明

### 图表导航

1. **图1 - 总体流程概览**: 从这里开始，了解整体架构
2. **图2 - 初始化阶段**: 详细了解网格和数据结构的创建
3. **图3 - 单粒子输运**: 核心模拟循环，包含边界和碰撞处理子流程
4. **图4 - 格林函数累积**: 传递函数的详细计算过程
5. **图5 - 裂变矩阵记录**: 稀疏矩阵的填充过程
6. **图6 - 最终化处理**: HDF5输出和归一化
7. **图7 - 批次管理**: 数据在批次间的合并机制
8. **图8 - 数据结构关系**: 理解各数据结构之间的关系
9. **图9 - 时序图**: 查看粒子与两个功能模块的交互时序
10. **图10 - 数据处理总览**: 从输入到输出的数据流

### 颜色编码

- 🟢 **绿色**: 开始节点
- 🔵 **蓝色**: 格林函数相关操作
- 🟣 **紫色**: 裂变矩阵相关操作
- 🟠 **橙色**: 原子操作/线程安全操作
- 🟡 **黄色**: 互斥锁操作
- 🩷 **粉色**: 结束节点

### 关键概念

- **source_particle_id**: 源粒子唯一标识，次级粒子继承该ID
- **稀疏矩阵**: 使用 `key = i × n + j` 的线性索引存储
- **COO格式**: 三数组存储 (row[], col[], data[])
- **原子操作**: 保证多线程安全的累积操作
- **批次管理**: 定期合并临时数据到永久存储

---

**文档版本**: v1.0  
**创建日期**: 2025年11月3日  
**适用于**: OpenMC传递函数与裂变矩阵功能
