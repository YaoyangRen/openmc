"""
可视化共轭通量分布
读取 adjoint_flux.h5 文件并绘制三维空间中的共轭通量分布
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors

def read_adjoint_flux(filename='adjoint_flux.h5'):
    """读取共轭通量 HDF5 文件"""
    
    print(f"读取文件: {filename}")
    
    with h5py.File(filename, 'r') as f:
        # 读取属性
        print("\n文件属性:")
        for key in f.attrs.keys():
            print(f"  {key}: {f.attrs[key]}")
        
        # 读取网格参数
        origin = f['origin'][:] # type: ignore
        shape = f['shape'][:] # pyright: ignore[reportIndexIssue]
        
        print(f"\n网格信息:")
        print(f"  原点: {origin}")
        print(f"  维度: {shape} (nx × ny × nz)")
        print(f"  总单元数: {np.prod(shape)}")
        
        # 读取稀疏格式数据
        cell_indices = f['cell_indices'][:] # type: ignore
        adjoint_flux_values = f['adjoint_flux_values'][:] # type: ignore
        
        print(f"\n稀疏数据:")
        print(f"  非零单元数: {len(cell_indices)}")
        print(f"  稀疏度: {(1 - len(cell_indices)/np.prod(shape))*100:.2f}%")
        
        # 读取稠密格式数据
        adjoint_flux_dense = f['adjoint_flux_dense'][:]
        
        # 重塑为 3D 数组
        nx, ny, nz = shape
        adjoint_flux_3d = adjoint_flux_dense.reshape(nx, ny, nz)
        
        print(f"\n共轭通量统计:")
        print(f"  最小值: {np.min(adjoint_flux_3d):.6e}")
        print(f"  最大值: {np.max(adjoint_flux_3d):.6e}")
        print(f"  平均值: {np.mean(adjoint_flux_3d[adjoint_flux_3d > 0]):.6e}")
        
        return {
            'origin': origin,
            'shape': shape,
            'pitch': f.attrs['pitch'] if 'pitch' in f.attrs else 1.0,
            'adjoint_flux_3d': adjoint_flux_3d,
            'cell_indices': cell_indices,
            'adjoint_flux_values': adjoint_flux_values
        }


def plot_2d_slices(data, output_file='adjoint_flux_slices.png'):
    """绘制 XY、XZ、YZ 三个方向的切片"""
    
    adjoint_flux_3d = data['adjoint_flux_3d']
    shape = data['shape']
    pitch = data['pitch']
    origin = data['origin']
    
    nx, ny, nz = shape
    
    # 选择中心切片
    ix_center = nx // 2
    iy_center = ny // 2
    iz_center = nz // 2
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 创建坐标轴
    x_coords = origin[0] + np.arange(nx) * pitch
    y_coords = origin[1] + np.arange(ny) * pitch
    z_coords = origin[2] + np.arange(nz) * pitch
    
    # XY 切片 (z = center)
    xy_slice = adjoint_flux_3d[:, :, iz_center]
    
    # 线性尺度
    im1 = axes[0, 0].imshow(
        xy_slice.T,
        origin='lower',
        cmap='hot',
        extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
        aspect='auto'
    )
    axes[0, 0].set_title(f'XY Slice (z = {z_coords[iz_center]:.2f} cm)')
    axes[0, 0].set_xlabel('x (cm)')
    axes[0, 0].set_ylabel('y (cm)')
    plt.colorbar(im1, ax=axes[0, 0], label='Adjoint Flux')
    
    # 对数尺度
    xy_log = np.log10(xy_slice + 1e-30)
    im2 = axes[1, 0].imshow(
        xy_log.T,
        origin='lower',
        cmap='viridis',
        extent=[x_coords[0], x_coords[-1], y_coords[0], y_coords[-1]],
        aspect='auto'
    )
    axes[1, 0].set_title(f'XY Slice - Log Scale')
    axes[1, 0].set_xlabel('x (cm)')
    axes[1, 0].set_ylabel('y (cm)')
    plt.colorbar(im2, ax=axes[1, 0], label='log₁₀(Φ†)')
    
    # XZ 切片 (y = center)
    xz_slice = adjoint_flux_3d[:, iy_center, :]
    
    im3 = axes[0, 1].imshow(
        xz_slice.T,
        origin='lower',
        cmap='hot',
        extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
        aspect='auto'
    )
    axes[0, 1].set_title(f'XZ Slice (y = {y_coords[iy_center]:.2f} cm)')
    axes[0, 1].set_xlabel('x (cm)')
    axes[0, 1].set_ylabel('z (cm)')
    plt.colorbar(im3, ax=axes[0, 1], label='Adjoint Flux')
    
    xz_log = np.log10(xz_slice + 1e-30)
    im4 = axes[1, 1].imshow(
        xz_log.T,
        origin='lower',
        cmap='viridis',
        extent=[x_coords[0], x_coords[-1], z_coords[0], z_coords[-1]],
        aspect='auto'
    )
    axes[1, 1].set_title(f'XZ Slice - Log Scale')
    axes[1, 1].set_xlabel('x (cm)')
    axes[1, 1].set_ylabel('z (cm)')
    plt.colorbar(im4, ax=axes[1, 1], label='log₁₀(Φ†)')
    
    # YZ 切片 (x = center)
    yz_slice = adjoint_flux_3d[ix_center, :, :]
    
    im5 = axes[0, 2].imshow(
        yz_slice.T,
        origin='lower',
        cmap='hot',
        extent=[y_coords[0], y_coords[-1], z_coords[0], z_coords[-1]],
        aspect='auto'
    )
    axes[0, 2].set_title(f'YZ Slice (x = {x_coords[ix_center]:.2f} cm)')
    axes[0, 2].set_xlabel('y (cm)')
    axes[0, 2].set_ylabel('z (cm)')
    plt.colorbar(im5, ax=axes[0, 2], label='Adjoint Flux')
    
    yz_log = np.log10(yz_slice + 1e-30)
    im6 = axes[1, 2].imshow(
        yz_log.T,
        origin='lower',
        cmap='viridis',
        extent=[y_coords[0], y_coords[-1], z_coords[0], z_coords[-1]],
        aspect='auto'
    )
    axes[1, 2].set_title(f'YZ Slice - Log Scale')
    axes[1, 2].set_xlabel('y (cm)')
    axes[1, 2].set_ylabel('z (cm)')
    plt.colorbar(im6, ax=axes[1, 2], label='log₁₀(Φ†)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n2D 切片图已保存: {output_file}")
    plt.close()


def plot_3d_scatter(data, output_file='adjoint_flux_3d_scatter.png', threshold_percentile=50):
    """绘制 3D 散点图（只显示高值区域）"""
    
    adjoint_flux_3d = data['adjoint_flux_3d']
    shape = data['shape']
    pitch = data['pitch']
    origin = data['origin']
    
    nx, ny, nz = shape
    
    # 创建坐标网格
    x = origin[0] + np.arange(nx) * pitch
    y = origin[1] + np.arange(ny) * pitch
    z = origin[2] + np.arange(nz) * pitch
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # 只绘制超过阈值的点
    threshold = np.percentile(adjoint_flux_3d[adjoint_flux_3d > 0], threshold_percentile)
    mask = adjoint_flux_3d > threshold
    
    print(f"\n3D 散点图:")
    print(f"  阈值 (>{threshold_percentile}%): {threshold:.6e}")
    print(f"  显示点数: {np.sum(mask)} / {np.prod(shape)}")
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制散点
    scatter = ax.scatter(
        X[mask], 
        Y[mask], 
        Z[mask],
        c=adjoint_flux_3d[mask],
        cmap='hot',
        s=20,
        alpha=0.6,
        norm=colors.LogNorm(vmin=threshold, vmax=np.max(adjoint_flux_3d))
    )
    
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    ax.set_zlabel('z (cm)')
    ax.set_title(f'Adjoint Flux 3D Distribution (>{threshold_percentile}th percentile)')
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('Adjoint Flux Φ†')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D 散点图已保存: {output_file}")
    plt.close()


def plot_3d_isosurface(data, output_file='adjoint_flux_3d_isosurface.png', num_levels=5):
    """绘制 3D 等值面投影"""
    
    adjoint_flux_3d = data['adjoint_flux_3d']
    shape = data['shape']
    pitch = data['pitch']
    origin = data['origin']
    
    nx, ny, nz = shape
    
    # 创建坐标
    x = origin[0] + np.arange(nx) * pitch
    y = origin[1] + np.arange(ny) * pitch
    z = origin[2] + np.arange(nz) * pitch
    
    fig = plt.figure(figsize=(18, 6))
    
    # 获取非零值的范围
    nonzero_values = adjoint_flux_3d[adjoint_flux_3d > 0]
    levels = np.logspace(
        np.log10(np.percentile(nonzero_values, 10)),
        np.log10(np.max(nonzero_values)),
        num_levels
    )
    
    # XY 投影 (最大值投影)
    ax1 = fig.add_subplot(131)
    xy_projection = np.max(adjoint_flux_3d, axis=2)
    
    X_xy, Y_xy = np.meshgrid(x, y, indexing='ij')
    contour1 = ax1.contourf(X_xy, Y_xy, xy_projection, levels=levels, cmap='hot', norm=colors.LogNorm())
    ax1.contour(X_xy, Y_xy, xy_projection, levels=levels, colors='black', alpha=0.3, linewidths=0.5)
    ax1.set_xlabel('x (cm)')
    ax1.set_ylabel('y (cm)')
    ax1.set_title('XY Projection (max along z)')
    ax1.set_aspect('equal')
    plt.colorbar(contour1, ax=ax1, label='Φ†')
    
    # XZ 投影
    ax2 = fig.add_subplot(132)
    xz_projection = np.max(adjoint_flux_3d, axis=1)
    
    X_xz, Z_xz = np.meshgrid(x, z, indexing='ij')
    contour2 = ax2.contourf(X_xz, Z_xz, xz_projection, levels=levels, cmap='hot', norm=colors.LogNorm())
    ax2.contour(X_xz, Z_xz, xz_projection, levels=levels, colors='black', alpha=0.3, linewidths=0.5)
    ax2.set_xlabel('x (cm)')
    ax2.set_ylabel('z (cm)')
    ax2.set_title('XZ Projection (max along y)')
    ax2.set_aspect('equal')
    plt.colorbar(contour2, ax=ax2, label='Φ†')
    
    # YZ 投影
    ax3 = fig.add_subplot(133)
    yz_projection = np.max(adjoint_flux_3d, axis=0)
    
    Y_yz, Z_yz = np.meshgrid(y, z, indexing='ij')
    contour3 = ax3.contourf(Y_yz, Z_yz, yz_projection, levels=levels, cmap='hot', norm=colors.LogNorm())
    ax3.contour(Y_yz, Z_yz, yz_projection, levels=levels, colors='black', alpha=0.3, linewidths=0.5)
    ax3.set_xlabel('y (cm)')
    ax3.set_ylabel('z (cm)')
    ax3.set_title('YZ Projection (max along x)')
    ax3.set_aspect('equal')
    plt.colorbar(contour3, ax=ax3, label='Φ†')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D 等值面投影已保存: {output_file}")
    plt.close()


def plot_sparse_distribution(data, output_file='adjoint_flux_sparse_dist.png'):
    """绘制稀疏数据的分布分析"""
    
    cell_indices = data['cell_indices']
    adjoint_flux_values = data['adjoint_flux_values']
    shape = data['shape']
    nx, ny, nz = shape
    
    # 将一维索引转换为 3D 坐标
    ix = cell_indices // (ny * nz)
    iy = (cell_indices % (ny * nz)) // nz
    iz = cell_indices % nz
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 值的直方图（对数尺度）
    axes[0, 0].hist(adjoint_flux_values, bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Adjoint Flux Value')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Distribution of Nonzero Values')
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 累积分布
    sorted_values = np.sort(adjoint_flux_values)[::-1]
    cumsum = np.cumsum(sorted_values)
    cumsum_normalized = cumsum / cumsum[-1]
    
    axes[0, 1].plot(np.arange(len(sorted_values)), cumsum_normalized, linewidth=2)
    axes[0, 1].axhline(y=0.9, color='r', linestyle='--', label='90%')
    axes[0, 1].axhline(y=0.95, color='orange', linestyle='--', label='95%')
    axes[0, 1].axhline(y=0.99, color='yellow', linestyle='--', label='99%')
    axes[0, 1].set_xlabel('Number of Cells (sorted by value)')
    axes[0, 1].set_ylabel('Cumulative Fraction')
    axes[0, 1].set_title('Cumulative Distribution of Adjoint Flux')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_xscale('log')
    
    # 3. X 方向分布
    axes[1, 0].hist(ix, bins=min(nx, 50), edgecolor='black', alpha=0.7, color='steelblue')
    axes[1, 0].set_xlabel('x index')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Spatial Distribution in X')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Y-Z 分布
    hist, xedges, yedges = np.histogram2d(iy, iz, bins=[min(ny, 30), min(nz, 30)])
    im = axes[1, 1].imshow(
        hist.T,
        origin='lower',
        cmap='YlOrRd',
        aspect='auto',
        extent=[0, ny, 0, nz]
    )
    axes[1, 1].set_xlabel('y index')
    axes[1, 1].set_ylabel('z index')
    axes[1, 1].set_title('Spatial Distribution in YZ Plane')
    plt.colorbar(im, ax=axes[1, 1], label='Count')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"稀疏数据分布分析已保存: {output_file}")
    plt.close()


def main():
    """主函数"""
    
    filename = 'adjoint_flux.h5'
    
    # 读取数据
    data = read_adjoint_flux(filename)
    
    # 生成各种可视化图
    print("\n开始生成可视化图像...")
    
    plot_2d_slices(data, 'adjoint_flux_slices.png')
    plot_3d_isosurface(data, 'adjoint_flux_3d_isosurface.png')
    plot_3d_scatter(data, 'adjoint_flux_3d_scatter.png', threshold_percentile=90)
    plot_sparse_distribution(data, 'adjoint_flux_sparse_dist.png')
    
    print("\n所有可视化完成！")
    print("生成的文件:")
    print("  - adjoint_flux_slices.png          (2D切片)")
    print("  - adjoint_flux_3d_isosurface.png   (3D等值面投影)")
    print("  - adjoint_flux_3d_scatter.png      (3D散点图)")
    print("  - adjoint_flux_sparse_dist.png     (稀疏数据分析)")


if __name__ == '__main__':
    main()
