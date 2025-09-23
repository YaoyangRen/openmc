#!/usr/bin/env python3
"""
详细分析格林函数HDF5文件的数据结构和存储格式
"""

import h5py
import numpy as np

def analyze_hdf5_structure():
    """详细分析HDF5文件结构"""
    print("=" * 70)
    print("格林函数HDF5文件结构分析")
    print("=" * 70)
    
    if not os.path.exists("green_function_data.h5"):
        print("错误: 未找到green_function_data.h5文件")
        return
    
    with h5py.File('green_function_data.h5', 'r') as f:
        print("\n1. 文件顶级结构:")
        print("-" * 40)
        
        def print_structure(name, obj):
            if isinstance(obj, h5py.Group):
                print(f"  📁 群组: {name}")
            elif isinstance(obj, h5py.Dataset):
                print(f"  📄 数据集: {name}")
                print(f"     类型: {obj.dtype}")
                print(f"     形状: {obj.shape}")
                print(f"     大小: {obj.size}")
        
        f.visititems(print_structure)
        
        print("\n2. 文件属性 (Attributes):")
        print("-" * 40)
        for attr_name in f.attrs.keys():
            attr_value = f.attrs[attr_name]
            if isinstance(attr_value, bytes):
                attr_value = attr_value.decode('utf-8')
            print(f"  {attr_name}: {attr_value}")
        
        print("\n3. 网格信息:")
        print("-" * 40)
        shape = f['shape'][...]
        origin = f['origin'][...]
        pitch = f.attrs['pitch']
        
        print(f"  网格形状 (shape): {shape}")
        print(f"  网格原点 (origin): {origin}")
        print(f"  网格间距 (pitch): {pitch}")
        
        # 计算网格边界
        end_point = origin + shape * pitch
        print(f"  网格终点: {end_point}")
        print(f"  网格覆盖范围:")
        for i, dim in enumerate(['X', 'Y', 'Z']):
            print(f"    {dim}: [{origin[i]:.2f}, {end_point[i]:.2f}] cm")
        
        total_cells = np.prod(shape)
        print(f"  总网格单元数: {total_cells}")
        
        print("\n4. 累积格林函数数据:")
        print("-" * 40)
        cumulative_gf = f['cumulative_green_function'][...]
        print(f"  数据类型: {cumulative_gf.dtype}")
        print(f"  数据形状: {cumulative_gf.shape}")
        print(f"  数据大小: {cumulative_gf.size}")
        print(f"  存储布局: 一维数组，按 index = ix + shape[0] * (iy + shape[1] * iz) 排列")
        print(f"  数值范围: [{np.min(cumulative_gf):.6e}, {np.max(cumulative_gf):.6e}]")
        print(f"  非零元素: {np.count_nonzero(cumulative_gf)} / {cumulative_gf.size}")
        print(f"  总和: {np.sum(cumulative_gf):.6e}")
        
        print("\n5. 源粒子数据:")
        print("-" * 40)
        particle_ids = f['source_particle_ids'][...]
        print(f"  源粒子ID数组:")
        print(f"    数据类型: {particle_ids.dtype}")
        print(f"    数组长度: {len(particle_ids)}")
        print(f"    ID范围: [{np.min(particle_ids)}, {np.max(particle_ids)}]")
        print(f"    前10个ID: {particle_ids[:10]}")
        
        print("\n6. 单个源粒子数据结构:")
        print("-" * 40)
        particles_group = f['source_particles']
        print(f"  源粒子群组中的数据集数量: {len(particles_group.keys())}")
        
        # 分析第一个粒子的数据
        first_particle_name = list(particles_group.keys())[0]
        first_particle_data = particles_group[first_particle_name][...]
        
        print(f"  示例粒子: {first_particle_name}")
        print(f"    数据类型: {first_particle_data.dtype}")
        print(f"    数据形状: {first_particle_data.shape}")
        print(f"    数据大小: {first_particle_data.size}")
        print(f"    数值范围: [{np.min(first_particle_data):.6e}, {np.max(first_particle_data):.6e}]")
        print(f"    非零元素: {np.count_nonzero(first_particle_data)} / {first_particle_data.size}")
        print(f"    总和: {np.sum(first_particle_data):.6e}")
        
        print("\n7. 数据访问示例:")
        print("-" * 40)
        print("  Python代码示例:")
        print("  ```python")
        print("  import h5py")
        print("  import numpy as np")
        print("  ")
        print("  with h5py.File('green_function_data.h5', 'r') as f:")
        print("      # 获取网格信息")
        print("      shape = f['shape'][...]")
        print("      origin = f['origin'][...]")
        print("      pitch = f.attrs['pitch']")
        print("      ")
        print("      # 获取累积格林函数")
        print("      cumulative_gf = f['cumulative_green_function'][...]")
        print("      gf_3d = cumulative_gf.reshape(shape)  # 转换为3D数组")
        print("      ")
        print("      # 获取源粒子ID列表")
        print("      particle_ids = f['source_particle_ids'][...]")
        print("      ")
        print("      # 获取特定源粒子的格林函数")
        print("      pid = particle_ids[0]")
        print("      particle_gf = f[f'source_particles/particle_{pid}'][...]")
        print("      particle_gf_3d = particle_gf.reshape(shape)")
        print("  ```")
        
        print("\n8. 存储效率分析:")
        print("-" * 40)
        file_size = os.path.getsize('green_function_data.h5')
        print(f"  文件总大小: {file_size / (1024**2):.2f} MB")
        
        # 估算理论大小
        n_particles = len(particle_ids)
        data_per_particle = cumulative_gf.size * 8  # 8 bytes per double
        total_particle_data = n_particles * data_per_particle
        cumulative_data_size = cumulative_gf.size * 8
        metadata_size = file_size - total_particle_data - cumulative_data_size
        
        print(f"  每个粒子数据大小: {data_per_particle / 1024:.2f} KB")
        print(f"  所有粒子数据大小: {total_particle_data / (1024**2):.2f} MB")
        print(f"  累积数据大小: {cumulative_data_size / 1024:.2f} KB")
        print(f"  元数据大小: {metadata_size / 1024:.2f} KB")
        
        print("\n9. 数据验证:")
        print("-" * 40)
        # 验证几个粒子的数据一致性
        total_from_particles = 0.0
        checked_particles = 0
        
        for i, pid in enumerate(particle_ids[:100]):  # 检查前100个粒子
            particle_data = particles_group[f'particle_{pid}'][...]
            total_from_particles += np.sum(particle_data)
            checked_particles += 1
        
        # 估算全部粒子的总和
        avg_contribution = total_from_particles / checked_particles
        estimated_total = avg_contribution * len(particle_ids)
        
        print(f"  检查的粒子数: {checked_particles}")
        print(f"  平均每粒子贡献: {avg_contribution:.6e}")
        print(f"  估算总贡献: {estimated_total:.6e}")
        print(f"  实际累积总和: {np.sum(cumulative_gf):.6e}")
        print(f"  相对误差: {abs(estimated_total - np.sum(cumulative_gf)) / np.sum(cumulative_gf) * 100:.4f}%")

if __name__ == "__main__":
    import os
    analyze_hdf5_structure()