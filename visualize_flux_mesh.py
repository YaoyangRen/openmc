#!/usr/bin/env python3
"""
通量网格可视化工具
Flux Mesh Visualization Tool

读取 flux_mesh.h5 文件并可视化空间通量分布
Reads flux_mesh.h5 and visualizes spatial flux distribution
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def read_flux_mesh(filename='flux_mesh.h5'):
    """
    读取通量网格 HDF5 文件
    Read flux mesh HDF5 file
    
    Parameters:
    -----------
    filename : str
        HDF5 文件路径
    
    Returns:
    --------
    dict : 包含网格数据和通量数据的字典
    """
    with h5py.File(filename, 'r') as f:
        # 读取网格参数
        grid_shape = f['grid_shape'][:]
        lower_left = f['grid_lower_left'][:]
        upper_right = f['grid_upper_right'][:]
        pitch = f['grid_pitch'][:]
        n_cells = f['n_cells'][()]
        n_batches = f['n_batches'][()]
        
        # 读取稀疏通量数据
        cell_indices = f['cell_indices'][:]
        flux_mean = f['flux_mean'][:]
        flux_std = f['flux_std'][:]
        
        # 读取密集格式数据
        flux_mean_dense = f['flux_mean_dense'][:]
        flux_std_dense = f['flux_std_dense'][:]
        
        data = {
            'grid_shape': grid_shape,
            'lower_left': lower_left,
            'upper_right': upper_right,
            'pitch': pitch,
            'n_cells': n_cells,
            'n_batches': n_batches,
            'cell_indices': cell_indices,
            'flux_mean': flux_mean,
            'flux_std': flux_std,
            'flux_mean_dense': flux_mean_dense,
            'flux_std_dense': flux_std_dense
        }
    
    return data

def reshape_flux_3d(data):
    """
    将一维通量数据重塑为三维数组
    Reshape 1D flux data to 3D array
    """
    nx, ny, nz = data['grid_shape']
    flux_3d = data['flux_mean_dense'].reshape((nz, ny, nx))
    error_3d = data['flux_std_dense'].reshape((nz, ny, nx))
    return flux_3d, error_3d

def plot_xy_slice(data, z_index=None, save_fig=True):
    """
    绘制 XY 平面切片
    Plot XY plane slice
    
    Parameters:
    -----------
    data : dict
        通量网格数据
    z_index : int
        Z 方向的索引 (默认为中间层)
    save_fig : bool
        是否保存图像
    """
    flux_3d, error_3d = reshape_flux_3d(data)
    nx, ny, nz = data['grid_shape']
    
    if z_index is None:
        z_index = nz // 2
    
    # 获取 XY 切片
    flux_slice = flux_3d[z_index, :, :]
    
    # 创建网格坐标
    x = np.linspace(data['lower_left'][0], data['upper_right'][0], nx)
    y = np.linspace(data['lower_left'][1], data['upper_right'][1], ny)
    X, Y = np.meshgrid(x, y)
    
    # 绘图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 只绘制非零通量
    flux_slice_masked = np.ma.masked_where(flux_slice == 0, flux_slice)
    
    im = ax.pcolormesh(X, Y, flux_slice_masked, cmap='jet', shading='auto')
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'Flux Distribution (XY Slice at z_index={z_index})', fontsize=14)
    ax.set_aspect('equal')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Flux (n/cm²/source)', fontsize=12)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(f'flux_xy_slice_z{z_index}.png', dpi=300)
        print(f'图像已保存: flux_xy_slice_z{z_index}.png')
    
    plt.show()

def plot_xz_slice(data, y_index=None, save_fig=True):
    """
    绘制 XZ 平面切片
    Plot XZ plane slice
    """
    flux_3d, error_3d = reshape_flux_3d(data)
    nx, ny, nz = data['grid_shape']
    
    if y_index is None:
        y_index = ny // 2
    
    # 获取 XZ 切片
    flux_slice = flux_3d[:, y_index, :]
    
    # 创建网格坐标
    x = np.linspace(data['lower_left'][0], data['upper_right'][0], nx)
    z = np.linspace(data['lower_left'][2], data['upper_right'][2], nz)
    X, Z = np.meshgrid(x, z)
    
    # 绘图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    flux_slice_masked = np.ma.masked_where(flux_slice == 0, flux_slice)
    
    im = ax.pcolormesh(X, Z, flux_slice_masked, cmap='jet', shading='auto')
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Z (cm)', fontsize=12)
    ax.set_title(f'Flux Distribution (XZ Slice at y_index={y_index})', fontsize=14)
    ax.set_aspect('equal')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Flux (n/cm²/source)', fontsize=12)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(f'flux_xz_slice_y{y_index}.png', dpi=300)
        print(f'图像已保存: flux_xz_slice_y{y_index}.png')
    
    plt.show()

def plot_relative_error(data, z_index=None, save_fig=True):
    """
    绘制相对误差分布
    Plot relative error distribution
    """
    flux_3d, error_3d = reshape_flux_3d(data)
    nx, ny, nz = data['grid_shape']
    
    if z_index is None:
        z_index = nz // 2
    
    # 计算相对误差
    flux_slice = flux_3d[z_index, :, :]
    error_slice = error_3d[z_index, :, :]
    
    # 避免除零
    relative_error = np.zeros_like(flux_slice)
    mask = flux_slice > 0
    relative_error[mask] = error_slice[mask] / flux_slice[mask] * 100
    
    # 创建网格坐标
    x = np.linspace(data['lower_left'][0], data['upper_right'][0], nx)
    y = np.linspace(data['lower_left'][1], data['upper_right'][1], ny)
    X, Y = np.meshgrid(x, y)
    
    # 绘图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    relative_error_masked = np.ma.masked_where(flux_slice == 0, relative_error)
    
    im = ax.pcolormesh(X, Y, relative_error_masked, cmap='viridis', shading='auto')
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'Relative Error (%) at z_index={z_index}', fontsize=14)
    ax.set_aspect('equal')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Relative Error (%)', fontsize=12)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(f'flux_error_z{z_index}.png', dpi=300)
        print(f'图像已保存: flux_error_z{z_index}.png')
    
    plt.show()

def plot_3d_scatter(data, threshold=None, save_fig=True):
    """
    绘制三维散点图 (仅显示非零通量的网格单元)
    Plot 3D scatter plot (only non-zero flux cells)
    
    Parameters:
    -----------
    data : dict
        通量网格数据
    threshold : float
        通量阈值 (仅显示大于此值的单元)
    save_fig : bool
        是否保存图像
    """
    flux_3d, _ = reshape_flux_3d(data)
    nx, ny, nz = data['grid_shape']
    
    # 创建三维坐标
    x = np.linspace(data['lower_left'][0], data['upper_right'][0], nx)
    y = np.linspace(data['lower_left'][1], data['upper_right'][1], ny)
    z = np.linspace(data['lower_left'][2], data['upper_right'][2], nz)
    
    # 找到所有非零单元
    indices = np.argwhere(flux_3d > 0)
    
    if threshold is not None:
        # 应用阈值
        threshold_mask = flux_3d[indices[:, 0], indices[:, 1], indices[:, 2]] > threshold
        indices = indices[threshold_mask]
    
    # 获取对应的坐标和通量值
    x_coords = x[indices[:, 2]]
    y_coords = y[indices[:, 1]]
    z_coords = z[indices[:, 0]]
    flux_values = flux_3d[indices[:, 0], indices[:, 1], indices[:, 2]]
    
    # 归一化通量值用于颜色映射
    norm = plt.Normalize(vmin=flux_values.min(), vmax=flux_values.max())
    colors = cm.jet(norm(flux_values))
    
    # 绘图
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    scatter = ax.scatter(x_coords, y_coords, z_coords, 
                        c=flux_values, cmap='jet', 
                        s=20, alpha=0.6)
    
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_zlabel('Z (cm)', fontsize=12)
    ax.set_title('3D Flux Distribution', fontsize=14)
    
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label('Flux (n/cm²/source)', fontsize=12)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig('flux_3d_scatter.png', dpi=300)
        print('图像已保存: flux_3d_scatter.png')
    
    plt.show()

def print_statistics(data):
    """
    打印通量统计信息
    Print flux statistics
    """
    flux_mean = data['flux_mean']
    flux_std = data['flux_std']
    
    print("\n" + "="*60)
    print("通量网格统计信息 / Flux Mesh Statistics")
    print("="*60)
    print(f"网格尺寸 / Grid Shape: {data['grid_shape']}")
    print(f"总单元数 / Total Cells: {data['n_cells']}")
    print(f"非零单元数 / Non-zero Cells: {len(flux_mean)}")
    print(f"批次数 / Number of Batches: {data['n_batches']}")
    print(f"网格间距 / Grid Pitch: {data['pitch']} cm")
    print(f"网格范围 / Grid Range:")
    print(f"  X: [{data['lower_left'][0]:.2f}, {data['upper_right'][0]:.2f}] cm")
    print(f"  Y: [{data['lower_left'][1]:.2f}, {data['upper_right'][1]:.2f}] cm")
    print(f"  Z: [{data['lower_left'][2]:.2f}, {data['upper_right'][2]:.2f}] cm")
    
    # 检查是否有通量数据
    if len(flux_mean) == 0:
        print("\n  警告 / WARNING:")
        print("  没有通量数据被统计!")
        print("  No flux data was tallied!")
        print("\n可能的原因 / Possible reasons:")
        print("  1. settings::flux_mesh_on 未启用")
        print("  2. 模拟未运行或未完成")
        print("  3. 粒子未经过网格区域")
        print("="*60 + "\n")
        return
    
    print(f"\n通量统计 / Flux Statistics:")
    print(f"  最大值 / Max: {flux_mean.max():.6e} n/cm²/source")
    print(f"  最小值 / Min: {flux_mean[flux_mean > 0].min():.6e} n/cm²/source")
    print(f"  平均值 / Mean: {flux_mean.mean():.6e} n/cm²/source")
    
    # 计算平均相对误差
    rel_error = np.zeros_like(flux_mean)
    mask = flux_mean > 0
    rel_error[mask] = flux_std[mask] / flux_mean[mask] * 100
    print(f"\n相对误差 / Relative Error:")
    print(f"  平均值 / Mean: {rel_error[mask].mean():.2f}%")
    print(f"  最大值 / Max: {rel_error[mask].max():.2f}%")
    print("="*60 + "\n")

def main():
    """
    主函数 - 完整可视化流程
    Main function - complete visualization workflow
    """
    import os
    
    # 检查文件是否存在
    filename = 'flux_mesh.h5'
    if not os.path.exists(filename):
        print(f"错误: 文件 {filename} 不存在!")
        print(f"Error: File {filename} does not exist!")
        return
    
    # 读取数据
    print(f"读取通量网格文件: {filename}")
    data = read_flux_mesh(filename)
    
    # 打印统计信息
    print_statistics(data)
    
    # 检查是否有数据可供可视化
    if len(data['flux_mean']) == 0:
        print("⚠️  无法生成可视化图形 - 没有通量数据")
        print("⚠️  Cannot generate plots - no flux data available")
        return
    
    # 绘制各种图形
    print("生成可视化图形...")
    
    # XY 切片 (中间层)
    plot_xy_slice(data, save_fig=True)
    
    # XZ 切片 (中间层)
    plot_xz_slice(data, save_fig=True)
    
    # 相对误差分布
    plot_relative_error(data, save_fig=True)
    
    # 3D 散点图
    # plot_3d_scatter(data, save_fig=True)  # 可选 - 数据量大时较慢
    
    print("\n可视化完成!")
    print("Visualization completed!")

if __name__ == '__main__':
    main()
