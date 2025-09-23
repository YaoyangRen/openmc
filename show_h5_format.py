#!/usr/bin/env python3
"""
使用h5dump风格的方式显示HDF5文件结构
"""

import h5py
import numpy as np

def show_h5_structure():
    """显示HDF5文件的详细结构"""
    print("=" * 80)
    print("格林函数数据文件 (green_function_data.h5) 结构详解")
    print("=" * 80)
    
    with h5py.File('green_function_data.h5', 'r') as f:
        print("文件根目录 (/)")
        print("├── 属性 (Attributes):")
        for name in sorted(f.attrs.keys()):
            value = f.attrs[name]
            if isinstance(value, bytes):
                value = value.decode()
            print(f"│   ├── {name}: {value}")
        
        print("├── 数据集 (Datasets):")
        datasets = [name for name in f.keys() if isinstance(f[name], h5py.Dataset)]
        for i, name in enumerate(sorted(datasets)):
            is_last = i == len(datasets) - 1
            prefix = "│   └──" if is_last else "│   ├──"
            dataset = f[name]
            print(f"{prefix} {name}")
            print(f"│       ├── 数据类型: {dataset.dtype}")
            print(f"│       ├── 形状: {dataset.shape}")
            print(f"│       └── 大小: {dataset.size} 元素")
        
        print("└── 群组 (Groups):")
        groups = [name for name in f.keys() if isinstance(f[name], h5py.Group)]
        for group_name in sorted(groups):
            group = f[group_name]
            print(f"    └── {group_name}/")
            
            # 显示群组中的数据集数量
            datasets_in_group = list(group.keys())
            print(f"        ├── 包含 {len(datasets_in_group)} 个数据集")
            
            # 显示几个示例
            for i, dataset_name in enumerate(datasets_in_group[:3]):
                dataset = group[dataset_name]
                prefix = "├──" if i < 2 else "└──"
                print(f"        {prefix} {dataset_name}")
                print(f"        │   ├── 类型: {dataset.dtype}")
                print(f"        │   └── 形状: {dataset.shape}")
            
            if len(datasets_in_group) > 3:
                print(f"        └── ... (还有 {len(datasets_in_group)-3} 个类似数据集)")

def explain_data_format():
    """详细说明数据格式"""
    print("\n" + "=" * 80)
    print("数据格式详细说明")
    print("=" * 80)
    
    with h5py.File('green_function_data.h5', 'r') as f:
        print("\n1. 文件头部信息")
        print("-" * 50)
        print("文件属性包含了格林函数网格的基本信息：")
        print("• filetype: 'green_function_mesh_per_particle' - 文件类型标识")
        print("• version: '1.0' - 数据格式版本")
        print(f"• pitch: {f.attrs['pitch']} - 网格单元的边长 (cm)")
        print(f"• n_source_particles: {f.attrs['n_source_particles']} - 源粒子总数")
        
        print("\n2. 网格几何信息")
        print("-" * 50)
        shape = f['shape'][...]
        origin = f['origin'][...]
        pitch = f.attrs['pitch']
        
        print(f"• shape: {shape} - 网格在X, Y, Z方向的单元数")
        print(f"• origin: {origin} - 网格原点坐标 (cm)")
        print(f"• 网格边界: X[{origin[0]}, {origin[0] + shape[0]*pitch}], "
              f"Y[{origin[1]}, {origin[1] + shape[1]*pitch}], "
              f"Z[{origin[2]}, {origin[2] + shape[2]*pitch}] cm")
        print(f"• 总网格单元数: {np.prod(shape)}")
        
        print("\n3. 数据存储布局")
        print("-" * 50)
        print("格林函数数据以一维数组形式存储，索引计算公式为：")
        print("  index = ix + shape[0] * (iy + shape[1] * iz)")
        print("其中：")
        print("  ix, iy, iz = 网格单元在X, Y, Z方向的索引 (从0开始)")
        print("  shape[0], shape[1], shape[2] = 各方向的网格数量")
        print()
        print("转换为3D数组：")
        print("  data_3d = data_1d.reshape(shape)")
        print("  value = data_3d[ix, iy, iz]")
        
        print("\n4. 累积格林函数 (cumulative_green_function)")
        print("-" * 50)
        cumulative = f['cumulative_green_function']
        print(f"• 数据类型: {cumulative.dtype} (双精度浮点数)")
        print(f"• 数组大小: {cumulative.size} (与总网格单元数相等)")
        print("• 物理意义: 所有源粒子贡献的叠加")
        print("• 单位: 取决于具体的物理量 (如通量、反应率等)")
        
        print("\n5. 源粒子ID列表 (source_particle_ids)")
        print("-" * 50)
        particle_ids = f['source_particle_ids']
        print(f"• 数据类型: {particle_ids.dtype} (64位整数)")
        print(f"• 数组长度: {particle_ids.size}")
        print("• 内容: 所有参与计算的源粒子的唯一ID")
        print("• 用途: 索引单个源粒子的格林函数数据")
        
        print("\n6. 单个源粒子数据 (source_particles/particle_XXXX)")
        print("-" * 50)
        particles_group = f['source_particles']
        first_particle = list(particles_group.keys())[0]
        first_data = particles_group[first_particle]
        
        print(f"• 存储位置: /source_particles/ 群组下")
        print(f"• 命名规则: particle_{{源粒子ID}}")
        print(f"• 数据类型: {first_data.dtype} (双精度浮点数)")
        print(f"• 数组大小: {first_data.size} (与累积数据相同)")
        print("• 物理意义: 单个源粒子对各网格单元的贡献")
        print("• 线性叠加: sum(所有单个粒子) ≈ 累积格林函数")
        
        print("\n7. 数据访问示例")
        print("-" * 50)
        print("Python代码示例：")
        print("""
import h5py
import numpy as np

with h5py.File('green_function_data.h5', 'r') as f:
    # 读取网格信息
    shape = f['shape'][...]
    origin = f['origin'][...]
    pitch = f.attrs['pitch']
    
    # 读取累积格林函数并转换为3D
    cumulative_1d = f['cumulative_green_function'][...]
    cumulative_3d = cumulative_1d.reshape(shape)
    
    # 读取源粒子ID列表
    particle_ids = f['source_particle_ids'][...]
    
    # 读取特定源粒子的数据
    pid = particle_ids[0]  # 选择第一个粒子
    particle_1d = f[f'source_particles/particle_{pid}'][...]
    particle_3d = particle_1d.reshape(shape)
    
    # 访问特定位置的值 (例如: ix=5, iy=5, iz=5)
    value = cumulative_3d[5, 5, 5]
    
    # 计算位置坐标
    x = origin[0] + 5 * pitch
    y = origin[1] + 5 * pitch  
    z = origin[2] + 5 * pitch
""")

if __name__ == "__main__":
    import os
    if os.path.exists('green_function_data.h5'):
        show_h5_structure()
        explain_data_format()
    else:
        print("错误: 未找到 green_function_data.h5 文件")