# 传递函数与裂变矩阵算法伪代码

**文档版本**: v1.0  
**日期**: 2025年10月30日  
**作者**: OpenMC 开发团队

---

## 目录

1. [传递函数算法](#1-传递函数算法)
2. [裂变矩阵算法](#2-裂变矩阵算法)
3. [数据结构定义](#3-数据结构定义)
4. [关键公式对照](#4-关键公式对照)

---

## 1. 传递函数算法

### 1.1 理论基础

**传递函数定义**：
$$T(P_0 \rightarrow r) = \int^{\infty}_0 \int^{4\pi}_0 \bar{\nu}\Sigma_f(r,E') \Phi(r,\Omega',E'|S_0)\,d\Omega'\,dE' \tag{1.1}$$

**物理意义**：从初始相空间点 $P_0$ 出发的源中子在空间位置 $r$ 处产生的平均裂变中子数。

---

### 1.2 初始化算法

```python
# ============================================================================
# 算法 1.1: 传递函数网格初始化
# ============================================================================
FUNCTION InitializeTransferFunctionMesh(pitch, max_batches, auto_bounds):
    INPUT:
        pitch          : 网格间距 (cm)
        max_batches    : 最大批次数
        auto_bounds    : 是否自动检测边界
    
    OUTPUT:
        mesh           : 传递函数网格对象
    
    BEGIN
        # 1. 确定空间边界
        IF auto_bounds THEN
            bbox ← GetRootUniverseBoundingBox()
            lower_bound ← [bbox.xmin, bbox.ymin, bbox.zmin]
            upper_bound ← [bbox.xmax, bbox.ymax, bbox.zmax]
            
            # 添加边距以避免边界效应
            margin ← pitch × 0.1
            FOR i ← 0 TO 2 DO
                lower_bound[i] ← lower_bound[i] - margin
                upper_bound[i] ← upper_bound[i] + margin
            END FOR
        ELSE
            lower_bound ← manual_lower_bound
            upper_bound ← manual_upper_bound
        END IF
        
        # 2. 计算网格维度
        FOR i ← 0 TO 2 DO
            delta[i] ← upper_bound[i] - lower_bound[i]
            shape[i] ← CEIL(delta[i] / pitch) + 1
        END FOR
        
        spatial_size ← shape[0] × shape[1] × shape[2]
        
        # 3. 初始化数据结构
        mesh.particle_transfer_functions ← EMPTY_MAP()     # 每个源粒子的传递函数
        mesh.current_batch_data ← EMPTY_MAP()              # 当前批次数据
        mesh.cumulative_data ← ARRAY(spatial_size, 0.0)    # 累积传递函数
        
        mesh.origin ← lower_bound
        mesh.shape ← shape
        mesh.pitch ← pitch
        mesh.inv_pitch ← 1.0 / pitch
        mesh.spatial_size ← spatial_size
        
        RETURN mesh
    END
END FUNCTION
```

---

### 1.3 路径积分累积算法

```python
# ============================================================================
# 算法 1.2: 粒子输运路径上的传递函数累积
# ============================================================================
FUNCTION AccumulateTransferFunctionOnPath(mesh, particle):
    INPUT:
        mesh           : 传递函数网格对象
        particle       : 当前粒子状态
    
    DESCRIPTION:
        在粒子输运路径上累积传递函数贡献
        contribution = w × Δs × ν̄Σf / Σt
    
    BEGIN
        # 1. 获取粒子属性
        position ← particle.position           # 当前位置 r
        weight ← particle.weight               # 粒子权重 w
        distance ← particle.distance           # 输运距离 Δs
        source_id ← particle.source_particle_id  # 源粒子ID
        
        # 2. 获取材料核截面
        IF particle.material == VOID THEN
            RETURN  # 真空中不累积
        END IF
        
        nu_sigma_f ← particle.macro_xs.nu_fission  # ν̄Σf
        sigma_t ← particle.macro_xs.total          # Σt
        
        # 3. 计算传递函数贡献
        IF sigma_t > 0 THEN
            # contribution = w × (Δs/Σt) × ν̄Σf
            #              = w × Δs × ν̄Σf / Σt
            contribution ← weight × distance × nu_sigma_f / sigma_t
        ELSE
            contribution ← 0.0
        END IF
        
        # 4. 只累积非零贡献
        IF contribution > 0 THEN
            # 4.1 计算空间网格索引
            index ← PositionToIndex(mesh, position)
            
            IF index >= 0 AND index < mesh.spatial_size THEN
                # 4.2 线程安全地获取或创建粒子数据
                LOCK(mesh.data_mutex)
                    IF source_id NOT IN mesh.current_batch_data THEN
                        mesh.current_batch_data[source_id] ← ARRAY(spatial_size, 0.0)
                    END IF
                    particle_data ← mesh.current_batch_data[source_id]
                UNLOCK(mesh.data_mutex)
                
                # 4.3 原子累积（线程安全）
                ATOMIC_ADD(particle_data[index], contribution)
                ATOMIC_ADD(mesh.cumulative_data[index], contribution)
            END IF
        END IF
    END
END FUNCTION
```

---

### 1.4 裂变事件累积算法

```python
# ============================================================================
# 算法 1.3: 裂变事件处的传递函数累积
# ============================================================================
FUNCTION AccumulateTransferFunctionAtFission(mesh, particle):
    INPUT:
        mesh           : 传递函数网格对象
        particle       : 发生裂变的粒子
    
    DESCRIPTION:
        在裂变点累积传递函数贡献
        contribution = w × ν̄Σf × σf / σt
    
    BEGIN
        # 1. 验证是裂变事件
        IF NOT particle.is_fission THEN
            RETURN
        END IF
        
        # 2. 获取核数据
        position ← particle.position
        weight ← particle.weight_last          # 碰撞前权重
        source_id ← particle.source_particle_id
        
        nuclide_xs ← particle.neutron_xs[particle.event_nuclide]
        nu_fission ← nuclide_xs.nu_fission     # ν̄
        sigma_f ← nuclide_xs.fission           # σf
        sigma_t ← nuclide_xs.total             # σt
        
        # 3. 计算裂变中子产生贡献
        IF sigma_t > 0 THEN
            # contribution = w × ν̄ × (σf/σt)
            contribution ← weight × nu_fission × sigma_f / sigma_t
        ELSE
            contribution ← 0.0
        END IF
        
        # 4. 累积到网格
        IF contribution > 0 THEN
            index ← PositionToIndex(mesh, position)
            
            IF index >= 0 AND index < mesh.spatial_size THEN
                LOCK(mesh.data_mutex)
                    IF source_id NOT IN mesh.current_batch_data THEN
                        mesh.current_batch_data[source_id] ← ARRAY(spatial_size, 0.0)
                    END IF
                    particle_data ← mesh.current_batch_data[source_id]
                UNLOCK(mesh.data_mutex)
                
                ATOMIC_ADD(particle_data[index], contribution)
                ATOMIC_ADD(mesh.cumulative_data[index], contribution)
            END IF
        END IF
    END
END FUNCTION
```

---

### 1.5 批次管理算法

```python
# ============================================================================
# 算法 1.4: 批次间数据合并
# ============================================================================
FUNCTION StartNewBatch(mesh, batch_id):
    INPUT:
        mesh           : 传递函数网格对象
        batch_id       : 新批次ID
    
    DESCRIPTION:
        保存当前批次数据并开始新批次
    
    BEGIN
        # 1. 保存上一批次的数据
        IF mesh.current_batch_id >= 0 THEN
            FOR EACH (source_id, batch_data) IN mesh.current_batch_data DO
                # 1.1 计算批次总贡献
                total_contribution ← SUM(batch_data)
                
                # 1.2 只保存非零数据
                IF total_contribution > 0 THEN
                    # 1.3 创建或获取该源粒子的累积传递函数
                    IF source_id NOT IN mesh.particle_transfer_functions THEN
                        mesh.particle_transfer_functions[source_id] ← 
                            ARRAY(mesh.spatial_size, 0.0)
                    END IF
                    
                    # 1.4 累积到该源粒子的总传递函数
                    particle_func ← mesh.particle_transfer_functions[source_id]
                    FOR i ← 0 TO mesh.spatial_size - 1 DO
                        particle_func[i] ← particle_func[i] + batch_data[i]
                    END FOR
                END IF
            END FOR
        END IF
        
        # 2. 开始新批次
        mesh.current_batch_id ← batch_id
        mesh.current_batch_data.CLEAR()
    END
END FUNCTION
```

---

### 1.6 终结与输出算法

```python
# ============================================================================
# 算法 1.5: 传递函数数据输出
# ============================================================================
FUNCTION FinalizeTransferFunction(mesh, filename):
    INPUT:
        mesh           : 传递函数网格对象
        filename       : 输出HDF5文件名
    
    OUTPUT:
        HDF5 文件包含所有源粒子的传递函数
    
    BEGIN
        # 1. 保存最后一个批次
        StartNewBatch(mesh, -1)
        
        # 2. 创建HDF5文件
        file ← CreateHDF5File(filename, WRITE_MODE)
        
        # 3. 写入元数据
        WriteAttribute(file, "filetype", "transfer_function_mesh_per_particle")
        WriteAttribute(file, "version", "2.0")
        WriteAttribute(file, "description", 
                      "Transfer function T(P0->r): Expected fission neutron production")
        WriteAttribute(file, "pitch", mesh.pitch)
        WriteDataset(file, "origin", mesh.origin)
        WriteDataset(file, "shape", mesh.shape)
        WriteAttribute(file, "n_source_particles", 
                      LENGTH(mesh.particle_transfer_functions))
        
        # 4. 写入累积传递函数（所有源粒子的总和）
        WriteDataset(file, "cumulative_transfer_function", mesh.cumulative_data)
        
        # 5. 为每个源粒子创建数据集
        particles_group ← CreateGroup(file, "source_particles")
        
        particle_ids ← EMPTY_ARRAY()
        FOR EACH (source_id, transfer_func) IN mesh.particle_transfer_functions DO
            # 5.1 写入该源粒子的传递函数
            dataset_name ← "particle_" + ToString(source_id)
            WriteDataset(particles_group, dataset_name, transfer_func)
            
            # 5.2 记录ID
            APPEND(particle_ids, source_id)
        END FOR
        
        # 6. 写入源粒子ID列表
        WriteDataset(file, "source_particle_ids", particle_ids)
        
        # 7. 关闭文件
        CloseGroup(particles_group)
        CloseFile(file)
        
        # 8. 输出统计信息
        PRINT("Transfer Function Finalization Complete")
        PRINT("  Grid shape: ", mesh.shape)
        PRINT("  Number of source particles: ", LENGTH(particle_ids))
        PRINT("  Output file: ", filename)
    END
END FUNCTION
```

---

## 2. 裂变矩阵算法

### 2.1 理论基础

**裂变矩阵定义**：
$$F[i][j] = \text{从空间单元 } i \text{ 出发的源中子在单元 } j \text{ 处产生的裂变中子数}$$

**数学表达**：
$$F[i][j] = \sum_{\text{源在}i} \frac{w \cdot \bar{\nu} \sigma_f}{\sigma_t \cdot k_{\text{eff}}}$$

---

### 2.2 初始化算法

```python
# ============================================================================
# 算法 2.1: 裂变矩阵初始化（稀疏存储）
# ============================================================================
FUNCTION InitializeFissionMatrix(pitch, max_batches, auto_bounds):
    INPUT:
        pitch          : 网格间距 (cm)
        max_batches    : 最大批次数
        auto_bounds    : 是否自动检测边界
    
    OUTPUT:
        matrix         : 裂变矩阵对象
    
    BEGIN
        # 1. 确定空间边界（同传递函数）
        IF auto_bounds THEN
            bbox ← GetRootUniverseBoundingBox()
            lower_bound ← [bbox.xmin - margin, bbox.ymin - margin, bbox.zmin - margin]
            upper_bound ← [bbox.xmax + margin, bbox.ymax + margin, bbox.zmax + margin]
        ELSE
            lower_bound ← manual_lower_bound
            upper_bound ← manual_upper_bound
        END IF
        
        # 2. 计算网格维度
        FOR i ← 0 TO 2 DO
            delta[i] ← upper_bound[i] - lower_bound[i]
            shape[i] ← CEIL(delta[i] / pitch) + 1
        END FOR
        
        n_cells ← shape[0] × shape[1] × shape[2]
        
        # 3. 初始化稀疏矩阵存储
        # 关键：使用哈希映射而非密集矩阵，节省95%+内存
        matrix.fission_matrix_sparse ← EMPTY_MAP()         # key=(i*n_cells+j), value=F[i][j]
        matrix.current_batch_sparse ← EMPTY_MAP()          # 当前批次的稀疏矩阵
        
        # 4. 源计数（用于归一化）
        matrix.source_counts ← ARRAY(n_cells, 0.0)         # 每个单元的源粒子数
        matrix.current_batch_source_counts ← ARRAY(n_cells, 0.0)
        
        # 5. 其他参数
        matrix.origin ← lower_bound
        matrix.shape ← shape
        matrix.pitch ← pitch
        matrix.inv_pitch ← 1.0 / pitch
        matrix.n_cells ← n_cells
        matrix.n_realizations ← 0
        
        RETURN matrix
    END
END FUNCTION
```

---

### 2.3 裂变事件记录算法

```python
# ============================================================================
# 算法 2.2: 裂变矩阵元素累积（稀疏存储）
# ============================================================================
FUNCTION RecordFissionEvent(matrix, source_position, fission_position, 
                            nu_fission_contribution):
    INPUT:
        matrix                   : 裂变矩阵对象
        source_position          : 源粒子出生位置
        fission_position         : 裂变发生位置
        nu_fission_contribution  : w × ν̄σf/σt 贡献值
    
    DESCRIPTION:
        记录一个裂变事件到稀疏矩阵
    
    BEGIN
        # 1. 计算源单元索引（行）
        i_source ← PositionToIndex(matrix, source_position)
        IF i_source < 0 THEN RETURN  # 源在边界外
        
        # 2. 计算裂变单元索引（列）
        j_fission ← PositionToIndex(matrix, fission_position)
        IF j_fission < 0 THEN RETURN  # 裂变在边界外
        
        # 3. 计算稀疏矩阵的线性键值
        # 使用行主序：key = row × n_cols + col
        key ← i_source × matrix.n_cells + j_fission
        
        # 4. 线程安全地累积到稀疏矩阵
        LOCK(matrix.data_mutex)
            # 4.1 累积裂变贡献（哈希映射自动处理不存在的键）
            matrix.current_batch_sparse[key] ← 
                matrix.current_batch_sparse.GET(key, 0.0) + nu_fission_contribution
            
            # 4.2 记录源计数（用于归一化）
            # 注意：每个源粒子只计数一次，这里可能需要标记
        UNLOCK(matrix.data_mutex)
    END
END FUNCTION
```

---

### 2.4 源粒子追踪算法

```python
# ============================================================================
# 算法 2.3: 源粒子出生位置记录
# ============================================================================
FUNCTION RecordSourceBirth(matrix, source_position):
    INPUT:
        matrix           : 裂变矩阵对象
        source_position  : 源粒子出生位置
    
    DESCRIPTION:
        记录源粒子出生位置（用于后续归一化）
    
    BEGIN
        # 1. 计算源单元索引
        i_source ← PositionToIndex(matrix, source_position)
        
        IF i_source >= 0 AND i_source < matrix.n_cells THEN
            # 2. 原子累积源计数
            ATOMIC_ADD(matrix.current_batch_source_counts[i_source], 1.0)
        END IF
    END
END FUNCTION
```

---

### 2.5 批次管理与归一化

```python
# ============================================================================
# 算法 2.4: 批次间裂变矩阵合并
# ============================================================================
FUNCTION StartNewBatchFissionMatrix(matrix, batch_id):
    INPUT:
        matrix         : 裂变矩阵对象
        batch_id       : 新批次ID
    
    DESCRIPTION:
        合并当前批次数据到累积裂变矩阵（稀疏存储）
    
    BEGIN
        # 1. 保存上一批次的数据
        IF matrix.current_batch_id >= 0 THEN
            # 1.1 合并稀疏矩阵（O(nnz) 复杂度，而非 O(n²)）
            FOR EACH (key, value) IN matrix.current_batch_sparse DO
                matrix.fission_matrix_sparse[key] ← 
                    matrix.fission_matrix_sparse.GET(key, 0.0) + value
            END FOR
            
            # 1.2 累积源计数
            FOR i ← 0 TO matrix.n_cells - 1 DO
                matrix.source_counts[i] ← matrix.source_counts[i] + 
                                           matrix.current_batch_source_counts[i]
            END FOR
            
            # 1.3 增加实现次数
            matrix.n_realizations ← matrix.n_realizations + 1
        END IF
        
        # 2. 开始新批次
        matrix.current_batch_id ← batch_id
        matrix.current_batch_sparse.CLEAR()
        matrix.current_batch_source_counts.FILL(0.0)
    END
END FUNCTION
```

---

### 2.6 稀疏矩阵输出算法

```python
# ============================================================================
# 算法 2.5: 裂变矩阵输出（COO格式）
# ============================================================================
FUNCTION FinalizeFissionMatrix(matrix, filename):
    INPUT:
        matrix         : 裂变矩阵对象
        filename       : 输出HDF5文件名
    
    OUTPUT:
        HDF5 文件包含稀疏裂变矩阵（COO格式）
    
    BEGIN
        # 1. 保存最后一个批次
        StartNewBatchFissionMatrix(matrix, -1)
        
        # 2. 转换稀疏矩阵为COO格式
        nnz ← LENGTH(matrix.fission_matrix_sparse)  # 非零元素数量
        
        row_indices ← ARRAY(nnz)      # 行索引
        col_indices ← ARRAY(nnz)      # 列索引
        data_raw ← ARRAY(nnz)         # 原始数据
        data_normalized ← ARRAY(nnz)  # 归一化数据
        
        index ← 0
        FOR EACH (key, value) IN matrix.fission_matrix_sparse DO
            # 2.1 从线性键恢复(i,j)索引
            i ← key / matrix.n_cells        # 行索引（源单元）
            j ← key % matrix.n_cells        # 列索引（裂变单元）
            
            row_indices[index] ← i
            col_indices[index] ← j
            data_raw[index] ← value
            
            # 2.2 归一化：F_norm[i][j] = F_raw[i][j] / (n_realizations × k_eff)
            # 这样 sum_j(F_norm[i][j]) 给出从单元i产生的下一代源粒子的期望分布
            IF matrix.n_realizations > 0 AND k_effective > 0 THEN
                data_normalized[index] ← value / (matrix.n_realizations × k_effective)
            ELSE
                data_normalized[index] ← value
            END IF
            
            index ← index + 1
        END FOR
        
        # 3. 创建HDF5文件
        file ← CreateHDF5File(filename, WRITE_MODE)
        
        # 4. 写入元数据
        WriteAttribute(file, "filetype", "fission_matrix_sparse")
        WriteAttribute(file, "version", "2.0")
        WriteAttribute(file, "storage_format", "COO (Coordinate)")
        WriteAttribute(file, "pitch", matrix.pitch)
        WriteDataset(file, "origin", matrix.origin)
        WriteDataset(file, "shape", matrix.shape)
        WriteAttribute(file, "n_cells", matrix.n_cells)
        WriteAttribute(file, "nnz", nnz)
        WriteAttribute(file, "n_realizations", matrix.n_realizations)
        WriteAttribute(file, "k_effective", k_effective)
        
        # 5. 计算稀疏度统计
        sparsity ← (1.0 - nnz / (matrix.n_cells × matrix.n_cells)) × 100.0
        WriteAttribute(file, "sparsity_percent", sparsity)
        
        # 6. 写入COO格式数据
        WriteDataset(file, "row_indices", row_indices)
        WriteDataset(file, "col_indices", col_indices)
        WriteDataset(file, "data_raw", data_raw)
        WriteDataset(file, "data_normalized", data_normalized)
        
        # 7. 写入源计数
        WriteDataset(file, "source_counts", matrix.source_counts)
        
        # 8. 内存统计
        dense_memory_MB ← (matrix.n_cells² × 8) / (1024² )  # 8 bytes per double
        sparse_memory_MB ← (nnz × 24) / (1024²)            # 24 bytes per entry (2 ints + 1 double)
        compression_ratio ← (1.0 - sparse_memory_MB / dense_memory_MB) × 100.0
        
        WriteAttribute(file, "dense_memory_MB", dense_memory_MB)
        WriteAttribute(file, "sparse_memory_MB", sparse_memory_MB)
        WriteAttribute(file, "memory_savings_percent", compression_ratio)
        
        # 9. 关闭文件
        CloseFile(file)
        
        # 10. 输出统计信息
        PRINT("Fission Matrix Finalization Complete (Sparse Storage)")
        PRINT("  Grid shape: ", matrix.shape)
        PRINT("  Total cells: ", matrix.n_cells)
        PRINT("  Non-zero elements: ", nnz)
        PRINT("  Sparsity: ", sparsity, "%")
        PRINT("  Memory savings: ", compression_ratio, "%")
        PRINT("  Dense matrix would be: ", dense_memory_MB, " MB")
        PRINT("  Sparse matrix is: ", sparse_memory_MB, " MB")
        PRINT("  Output file: ", filename)
    END
END FUNCTION
```

---

## 3. 数据结构定义

### 3.1 传递函数网格结构

```python
STRUCTURE TransferFunctionMesh:
    # 网格几何
    origin : Array[3, double]              # 网格原点 [x0, y0, z0]
    shape : Array[3, int]                  # 网格形状 [nx, ny, nz]
    pitch : double                         # 网格间距
    inv_pitch : double                     # 间距倒数（优化计算）
    spatial_size : int                     # 总单元数 = nx × ny × nz
    
    # 传递函数数据（分粒子存储）
    particle_transfer_functions : Map<int64, Array[spatial_size, double>>
        # key: source_particle_id
        # value: 该源粒子的传递函数 T_i(P_i -> r)
    
    # 当前批次数据
    current_batch_data : Map<int64, Array[spatial_size, double>>
        # key: source_particle_id
        # value: 当前批次该粒子的累积
    
    # 累积所有源粒子的传递函数
    cumulative_data : Array[spatial_size, double]
    
    # 批次管理
    current_batch_id : int
    max_batches : int
    
    # 线程同步
    data_mutex : Mutex
    
    # 统计信息
    total_contributions : atomic<uint64>
    dropped_contributions : atomic<uint64>
END STRUCTURE
```

---

### 3.2 裂变矩阵结构（稀疏存储）

```python
STRUCTURE FissionMatrix:
    # 网格几何
    origin : Array[3, double]              # 网格原点
    shape : Array[3, int]                  # 网格形状 [nx, ny, nz]
    pitch : double                         # 网格间距
    inv_pitch : double                     # 间距倒数
    n_cells : int                          # 总单元数 = nx × ny × nz
    
    # 稀疏矩阵存储（COO格式的哈希表实现）
    fission_matrix_sparse : Map<size_t, double>
        # key = row × n_cells + col  (线性索引)
        # value = F[row][col] (裂变矩阵元素)
        # 只存储非零元素，节省95-99%内存
    
    # 当前批次稀疏矩阵
    current_batch_sparse : Map<size_t, double>
    
    # 源粒子计数（用于归一化）
    source_counts : Array[n_cells, double]
        # source_counts[i] = 从单元i出发的源粒子总数
    
    current_batch_source_counts : Array[n_cells, double]
    
    # 统计信息
    n_realizations : int                   # 实现次数（批次数）
    
    # 批次管理
    current_batch_id : int
    max_batches : int
    
    # 线程同步
    data_mutex : Mutex
END STRUCTURE
```

---

### 3.3 粒子状态结构

```python
STRUCTURE Particle:
    # 位置和方向
    position : Position                    # 当前位置 (x, y, z)
    direction : Direction                  # 运动方向 (u, v, w)
    
    # 能量和权重
    energy : double                        # 中子能量 (MeV)
    weight : double                        # 当前权重
    weight_last : double                   # 碰撞前权重
    
    # 输运信息
    distance : double                      # 本次输运距离
    material : int                         # 当前材料ID
    
    # 核截面数据
    macro_xs : MacroscopicXS
        total : double                     # 总截面 Σt
        fission : double                   # 裂变截面 Σf
        nu_fission : double                # ν̄Σf
        absorption : double                # 吸收截面
        scattering : double                # 散射截面
    
    # 事件信息
    event_type : EventType                 # SCATTER, FISSION, ABSORPTION, etc.
    event_nuclide : int                    # 发生反应的核素索引
    is_fission : bool                      # 是否发生裂变
    
    # 源粒子追踪
    source_particle_id : int64             # 所属源粒子的唯一ID
    source_position : Position             # 源粒子出生位置
    
    # 类型
    type : ParticleType                    # NEUTRON, PHOTON, etc.
END STRUCTURE
```

---

## 4. 关键公式对照

### 4.1 传递函数公式对照表

| 理论公式 | 蒙特卡洛实现 | 说明 |
|---------|------------|------|
| $T(P_0 \rightarrow r) = \int\int \bar{\nu}\Sigma_f \Phi \,d\Omega dE$ | `contribution = w × Δs × ν̄Σf / Σt` | 路径积分累积 |
| $\Phi(r,E)$ (通量) | `w × Δs / Σt` | track-length estimator |
| $\bar{\nu}\Sigma_f(r,E) \Phi(r,E)$ | `(w × Δs / Σt) × ν̄Σf` | 裂变中子产生率 |
| 裂变点贡献 | `contribution = w × ν̄ × σf / σt` | 事件发生时 |

---

### 4.2 裂变矩阵公式对照表

| 理论定义 | 蒙特卡洛实现 | 说明 |
|---------|------------|------|
| $F[i][j]$ | `fission_matrix_sparse[i × n + j]` | 稀疏存储 |
| 从单元$i$到单元$j$的裂变 | `i = src_index, j = fis_index` | 空间映射 |
| 裂变中子贡献 | `w × ν̄σf / σt` | 单次事件 |
| 归一化 | `F[i][j] / (n_batches × k_eff)` | 期望估计 |
| 稀疏度 | `1 - nnz / n²` | 通常 > 95% |

---

### 4.3 空间索引映射

```python
# ============================================================================
# 算法 4.1: 位置到网格索引的转换
# ============================================================================
FUNCTION PositionToIndex(mesh, position):
    INPUT:
        mesh       : 网格对象（传递函数或裂变矩阵）
        position   : 三维位置 (x, y, z)
    
    OUTPUT:
        index      : 线性索引（0 到 spatial_size - 1），或 -1 表示超出边界
    
    BEGIN
        eps ← 1.0e-10  # 数值容差
        
        # 1. 计算三维网格索引
        ix ← FLOOR((position.x - mesh.origin[0] + eps) × mesh.inv_pitch)
        iy ← FLOOR((position.y - mesh.origin[1] + eps) × mesh.inv_pitch)
        iz ← FLOOR((position.z - mesh.origin[2] + eps) × mesh.inv_pitch)
        
        # 2. 边界检查
        IF ix < 0 OR ix >= mesh.shape[0] THEN RETURN -1
        IF iy < 0 OR iy >= mesh.shape[1] THEN RETURN -1
        IF iz < 0 OR iz >= mesh.shape[2] THEN RETURN -1
        
        # 3. 转换为线性索引（行主序）
        # index = ix + nx × (iy + ny × iz)
        index ← ix + mesh.shape[0] × (iy + mesh.shape[1] × iz)
        
        RETURN index
    END
END FUNCTION
```

---

## 5. 完整工作流程

### 5.1 传递函数计算流程

```
初始化阶段:
    ├─ InitializeTransferFunctionMesh()
    │   ├─ 确定空间边界
    │   ├─ 创建网格
    │   └─ 初始化数据结构
    │
模拟阶段 (每个批次):
    ├─ 对每个源粒子:
    │   ├─ 记录 source_particle_id 和出生位置
    │   ├─ 输运过程:
    │   │   ├─ 每次移动调用 AccumulateTransferFunctionOnPath()
    │   │   │   └─ 累积 w × Δs × ν̄Σf / Σt
    │   │   └─ 裂变事件调用 AccumulateTransferFunctionAtFission()
    │   │       └─ 累积 w × ν̄ × σf / σt
    │   └─ 粒子结束
    │
    └─ 批次结束调用 StartNewBatch()
        └─ 合并批次数据到 particle_transfer_functions
        
终结阶段:
    └─ FinalizeTransferFunction()
        ├─ 保存最后批次数据
        ├─ 转换为HDF5格式
        ├─ 写入每个源粒子的传递函数
        └─ 写入累积传递函数
```

---

### 5.2 裂变矩阵计算流程

```
初始化阶段:
    ├─ InitializeFissionMatrix()
    │   ├─ 确定空间边界
    │   ├─ 创建网格
    │   └─ 初始化稀疏矩阵（哈希映射）
    │
模拟阶段 (每个活跃批次):
    ├─ 对每个源粒子:
    │   ├─ RecordSourceBirth(源位置)
    │   │   └─ source_counts[i_source] += 1
    │   │
    │   ├─ 输运过程中，每次裂变:
    │   │   └─ RecordFissionEvent(源位置, 裂变位置, ν̄Σf贡献)
    │   │       ├─ i = PositionToIndex(源位置)
    │   │       ├─ j = PositionToIndex(裂变位置)
    │   │       ├─ key = i × n_cells + j
    │   │       └─ current_batch_sparse[key] += contribution
    │   │
    │   └─ 粒子结束
    │
    └─ 批次结束调用 StartNewBatchFissionMatrix()
        ├─ 合并 current_batch_sparse 到 fission_matrix_sparse
        └─ 累积 source_counts
        
终结阶段:
    └─ FinalizeFissionMatrix()
        ├─ 保存最后批次数据
        ├─ 转换稀疏矩阵为COO格式
        │   ├─ 提取 (row, col, value) 三元组
        │   └─ 归一化: value / (n_batches × k_eff)
        ├─ 写入HDF5文件
        │   ├─ row_indices[], col_indices[]
        │   ├─ data_raw[], data_normalized[]
        │   └─ 元数据和统计信息
        └─ 输出内存节省统计
```

---

## 6. 性能优化要点

### 6.1 传递函数优化

| 优化技术 | 实现方式 | 效果 |
|---------|---------|------|
| 线程安全 | `std::mutex` + 原子操作 | 支持OpenMP并行 |
| 早期退出 | `if (contribution > 0)` | 减少无效写入 |
| 内存预分配 | 批次开始时分配 | 减少动态分配 |
| 分粒子存储 | `Map<id, data>` | 精确追踪每个源粒子 |

---

### 6.2 裂变矩阵优化

| 优化技术 | 实现方式 | 效果 |
|---------|---------|------|
| 稀疏存储 | `Map<key, value>` | 内存节省 95-99% |
| 线性索引 | `key = i × n + j` | O(1) 查找 |
| 批次合并 | 延迟合并到全局 | 减少锁竞争 |
| COO输出 | 三数组格式 | 标准稀疏格式 |

---

### 6.3 复杂度分析

**传递函数**:

- 空间复杂度: $O(N_p \times N_{\text{cells}})$，其中 $N_p$ 是源粒子数
- 时间复杂度: $O(N_{\text{tracks}})$ 每次累积 $O(1)$

**裂变矩阵（稀疏）**:

- 空间复杂度: $O(\text{nnz})$，其中 $\text{nnz} \ll N_{\text{cells}}^2$
- 时间复杂度: $O(\text{nnz})$ 批次合并，而非 $O(N_{\text{cells}}^2)$

**内存对比**:

```
密集裂变矩阵: n² × 8 bytes
稀疏裂变矩阵: nnz × 24 bytes
节省比例: 1 - (nnz × 24) / (n² × 8)

例子 (50³ = 125,000 cells, sparsity = 99%):
  密集: 125,000² × 8 ≈ 117 GB
  稀疏: 125,000 × 24 ≈ 35.7 MB
  节省: 99.97%
```

---

## 7. 参考文献

1. **传递函数理论**:
   - Equation (3.2.1): $T(P_0 \rightarrow r) = \int \bar{\nu}\Sigma_f \Phi \,d\Omega dE$

2. **蒙特卡洛方法**:
   - Track-length flux estimator: $\Phi = w \times \Delta s / \Sigma_t$
   - Collision estimator: $\Phi = w / \Sigma_t$

3. **稀疏矩阵存储**:
   - COO (Coordinate) format: `(row_indices[], col_indices[], data[])`
   - 与 SciPy, MATLAB, Eigen 兼容

4. **裂变矩阵应用**:
   - 增殖因子计算: $k_{\text{eff}}$ = 最大特征值
   - 伴随通量求解
   - 扰动理论分析

---

**文档结束**

如需更多细节，请参考：

- `src/greenfunction_mesh.cpp` - 传递函数实现
- `src/fission_matrix.cpp` - 裂变矩阵实现
- `markdown/SPARSE_FISSION_MATRIX.md` - 稀疏矩阵详细文档
