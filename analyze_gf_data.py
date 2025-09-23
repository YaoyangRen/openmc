#!/usr/bin/env python3
"""
直接分析现有的格林函数数据文件
"""

import h5py
import numpy as np

def analyze_existing_data():
    """分析现有的格林函数数据"""
    print("分析现有的格林函数数据...")
    
    if not os.path.exists("green_function_data.h5"):
        print("错误: 未找到green_function_data.h5文件")
        return False
    
    try:
        with h5py.File('green_function_data.h5', 'r') as f:
            print(f"文件类型: {f.attrs['filetype']}")
            print(f"版本: {f.attrs['version']}")
            print(f"网格间距: {f.attrs['pitch']}")
            print(f"源粒子数量: {f.attrs['n_source_particles']}")
            
            # 获取网格信息
            shape = f['shape'][:]
            origin = f['origin'][:]
            print(f"网格形状: {shape}")
            print(f"网格原点: {origin}")
            
            # 获取源粒子ID列表
            particle_ids = f['source_particle_ids'][:]
            print(f"源粒子ID数量: {len(particle_ids)}")
            print(f"前10个源粒子ID: {particle_ids[:10]}")
            
            # 获取累积格林函数
            cumulative_gf = f['cumulative_green_function'][:]
            print(f"累积格林函数总和: {np.sum(cumulative_gf)}")
            print(f"累积格林函数最大值: {np.max(cumulative_gf)}")
            print(f"累积格林函数非零元素数: {np.count_nonzero(cumulative_gf)}")
            
            # 分析每个源粒子的贡献
            total_contribution = 0.0
            particle_contributions = []
            
            for pid in particle_ids:
                particle_name = f'source_particles/particle_{pid}'
                if particle_name in f:
                    particle_gf = f[particle_name][:]
                    contribution = np.sum(particle_gf)
                    particle_contributions.append((pid, contribution))
                    total_contribution += contribution
                else:
                    print(f"警告: 未找到粒子 {pid} 的数据")
            
            print(f"所有源粒子贡献总和: {total_contribution}")
            
            # 验证线性叠加性质
            diff = abs(total_contribution - np.sum(cumulative_gf))
            print(f"累积格林函数与单个贡献之差: {diff}")
            
            if diff < 1e-10:
                print("✓ 格林函数线性叠加性质验证通过!")
            else:
                print(f"✗ 格林函数线性叠加性质验证失败! 差值: {diff}")
            
            # 显示前几个源粒子的贡献
            sorted_contributions = sorted(particle_contributions, 
                                        key=lambda x: x[1], reverse=True)
            print("\n前5个源粒子的贡献:")
            for i, (pid, contrib) in enumerate(sorted_contributions[:5]):
                percentage = contrib / total_contribution * 100 if total_contribution > 0 else 0
                print(f"  源粒子 {pid}: {contrib:.6e} ({percentage:.2f}%)")
        
        return True
        
    except Exception as e:
        print(f"分析格林函数数据时出错: {e}")
        return False

if __name__ == "__main__":
    import os
    analyze_existing_data()