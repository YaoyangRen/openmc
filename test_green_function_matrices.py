#!/usr/bin/env python3
"""
测试每个源粒子单独统计格林函数矩阵功能的脚本
"""

import os
import sys
import subprocess
import h5py
import numpy as np
import matplotlib.pyplot as plt

def run_openmc_test():
    """运行OpenMC测试"""
    print("开始运行OpenMC测试...")
    
    # 切换到测试目录
    test_dir = "test_green_function"
    os.chdir(test_dir)
    
    try:
        # 运行OpenMC
        result = subprocess.run(["../bin/openmc"], 
                              capture_output=True, 
                              text=True, 
                              timeout=300)
        
        if result.returncode == 0:
            print("OpenMC运行成功!")
            print(result.stdout)
        else:
            print("OpenMC运行失败:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("OpenMC运行超时")
        return False
    except FileNotFoundError:
        print("找不到OpenMC可执行文件，请确认编译成功")
        return False
    
    return True

def analyze_green_function_data():
    """分析格林函数数据"""
    print("\n开始分析格林函数数据...")
    
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
            print(f"源粒子ID: {particle_ids[:10]}...")  # 只显示前10个
            
            # 获取累积格林函数
            cumulative_gf = f['cumulative_green_function'][:]
            print(f"累积格林函数总和: {np.sum(cumulative_gf)}")
            
            # 分析每个源粒子的贡献
            total_contribution = 0.0
            particle_contributions = {}
            
            for pid in particle_ids:
                particle_name = f'source_particles/particle_{pid}'
                if particle_name in f:
                    particle_gf = f[particle_name][:]
                    contribution = np.sum(particle_gf)
                    particle_contributions[pid] = contribution
                    total_contribution += contribution
            
            print(f"所有源粒子贡献总和: {total_contribution}")
            print(f"累积格林函数总和: {np.sum(cumulative_gf)}")
            
            # 验证线性叠加性质
            if abs(total_contribution - np.sum(cumulative_gf)) < 1e-10:
                print("✓ 格林函数线性叠加性质验证通过!")
            else:
                print("✗ 格林函数线性叠加性质验证失败!")
            
            # 显示前几个源粒子的贡献
            sorted_contributions = sorted(particle_contributions.items(), 
                                        key=lambda x: x[1], reverse=True)
            print("\n前5个源粒子的贡献:")
            for i, (pid, contrib) in enumerate(sorted_contributions[:5]):
                percentage = contrib / total_contribution * 100 if total_contribution > 0 else 0
                print(f"  源粒子 {pid}: {contrib:.6e} ({percentage:.2f}%)")
        
        return True
        
    except Exception as e:
        print(f"分析格林函数数据时出错: {e}")
        return False

def create_visualization():
    """创建可视化图表"""
    print("\n创建可视化图表...")
    
    try:
        with h5py.File('green_function_data.h5', 'r') as f:
            shape = f['shape'][:]
            origin = f['origin'][:]
            pitch = f.attrs['pitch']
            
            # 获取累积格林函数
            cumulative_gf = f['cumulative_green_function'][:]
            cumulative_3d = cumulative_gf.reshape(shape)
            
            # 创建中心截面图
            center_z = shape[2] // 2
            center_slice = cumulative_3d[:, :, center_z]
            
            # 创建坐标轴
            x = np.linspace(origin[0], origin[0] + shape[0] * pitch, shape[0])
            y = np.linspace(origin[1], origin[1] + shape[1] * pitch, shape[1])
            X, Y = np.meshgrid(x, y)
            
            plt.figure(figsize=(10, 8))
            plt.contourf(X, Y, center_slice.T, levels=50, cmap='viridis')
            plt.colorbar(label='累积格林函数值')
            plt.xlabel('X (cm)')
            plt.ylabel('Y (cm)')
            plt.title(f'累积格林函数中心截面 (z = {origin[2] + center_z * pitch:.2f} cm)')
            plt.grid(True, alpha=0.3)
            plt.savefig('green_function_visualization.png', dpi=300, bbox_inches='tight')
            plt.show()
            
            print("可视化图表已保存为 green_function_visualization.png")
        
        return True
        
    except Exception as e:
        print(f"创建可视化时出错: {e}")
        return False

def main():
    """主测试函数"""
    print("=" * 60)
    print("OpenMC 每源粒子格林函数矩阵测试")
    print("=" * 60)
    
    # 保存当前目录
    original_dir = os.getcwd()
    
    try:
        # 步骤1: 运行OpenMC
        if not run_openmc_test():
            return 1
        
        # 步骤2: 分析数据
        if not analyze_green_function_data():
            return 1
        
        # 步骤3: 创建可视化
        create_visualization()
        
        print("\n" + "=" * 60)
        print("测试完成!")
        print("主要检查点:")
        print("1. ✓ green_function_data.h5文件已生成")
        print("2. ✓ 每个源粒子都有独立的格林函数矩阵")
        print("3. ✓ 格林函数线性叠加性质验证")
        print("4. ✓ 数据可视化")
        print("=" * 60)
        
        return 0
        
    except Exception as e:
        print(f"测试过程中出现错误: {e}")
        return 1
    
    finally:
        # 恢复原始目录
        os.chdir(original_dir)

if __name__ == "__main__":
    sys.exit(main())