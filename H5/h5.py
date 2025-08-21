import h5py
import numpy as np

with h5py.File('green_function_data.h5', 'r') as f:
    print("文件属性:")
    print(f"  文件类型: {f.attrs['filetype']}")
    print(f"  网格间距: {f.attrs['pitch']}")
    print(f"  batch数量: {f.attrs['n_batches']}")
    
    # 读取网格参数
    origin = f['origin'][:]
    shape = f['shape'][:]
    print(f"  网格形状: {shape}")
    print(f"  网格原点: {origin}")
    
    # 分析每个batch的数据
    print(f"\n每代（Batch）统计数据:")
    print("-" * 80)
    print(f"{'Batch':>5} {'总和':>15} {'非零点数':>10} {'最大值':>15} {'平均值':>15}")
    print("-" * 80)
    
    batch_stats = []
    for i in range(f.attrs['n_batches']):
        batch_data = f[f'batch_{i}'][:]
        total_sum = np.sum(batch_data)
        nonzero_count = np.count_nonzero(batch_data)
        max_value = np.max(batch_data)
        mean_value = np.mean(batch_data[batch_data > 0]) if nonzero_count > 0 else 0
        
        print(f"{i:>5} {total_sum:>15.2e} {nonzero_count:>10} {max_value:>15.2e} {mean_value:>15.2e}")
        
        batch_stats.append({
            'batch': i,
            'total': total_sum,
            'nonzero_count': nonzero_count,
            'max_value': max_value,
            'mean_value': mean_value
        })
    
    print("-" * 80)
    
    # 总体统计
    batch_totals = [stat['total'] for stat in batch_stats]
    print(f"\n总体统计:")
    print(f"  总batch数: {len(batch_stats)}")
    print(f"  平均总值: {np.mean(batch_totals):.2e}")
    print(f"  标准差: {np.std(batch_totals):.2e}")
    print(f"  变异系数: {np.std(batch_totals)/np.mean(batch_totals)*100:.2f}%")
    print(f"  最小值: {np.min(batch_totals):.2e} (Batch {np.argmin(batch_totals)})")
    print(f"  最大值: {np.max(batch_totals):.2e} (Batch {np.argmax(batch_totals)})")
    
    # 可视化batch间变化
    plt.figure(figsize=(12, 8))
    
    # 第一个子图：batch变化
    plt.subplot(2, 2, 1)
    plt.plot(batch_totals)
    plt.xlabel('Batch Number')
    plt.ylabel('Total Green Function Value')
    plt.title('Green Function Value vs Batch')
    plt.grid(True)
    
    # 读取累积数据进行3D可视化
    cumulative = f['cumulative_data'][:]
    shape = f['shape'][:]
    print(f"\n网格数据形状: {cumulative.shape}")
    
    # 显示数据表格（前10x10的切片）
    print(f"\n累积数据前10x10x1切片:")
    if len(cumulative.shape) == 3:
        slice_data = cumulative[:min(10, shape[0]), :min(10, shape[1]), 0]
        print("X\\Y", end="")
        for j in range(slice_data.shape[1]):
            print(f"{j:>12}", end="")
        print()
        for i in range(slice_data.shape[0]):
            print(f"{i:>3}", end="")
            for j in range(slice_data.shape[1]):
                print(f"{slice_data[i,j]:>12.2e}", end="")
            print()
    
    # 2D热图 (Z=0切片)
    plt.subplot(2, 2, 2)
    if len(cumulative.shape) == 3:
        z_slice = cumulative[:, :, shape[2]//2]  # 中间Z切片
        im = plt.imshow(z_slice, origin='lower', aspect='auto', cmap='viridis')
        plt.colorbar(im)
        plt.title(f'Z={shape[2]//2} 切片热图')
        plt.xlabel('Y方向')
        plt.ylabel('X方向')
    
    # 3D散点图（只显示非零点）
    plt.subplot(2, 2, 3)
    if len(cumulative.shape) == 3:
        # 找到非零点
        nonzero_indices = np.nonzero(cumulative)
        if len(nonzero_indices[0]) > 0:
            # 限制显示点数以避免过密
            max_points = 1000
            if len(nonzero_indices[0]) > max_points:
                idx = np.random.choice(len(nonzero_indices[0]), max_points, replace=False)
                x_coords = nonzero_indices[0][idx]
                y_coords = nonzero_indices[1][idx]
                z_coords = nonzero_indices[2][idx]
                values = cumulative[x_coords, y_coords, z_coords]
            else:
                x_coords = nonzero_indices[0]
                y_coords = nonzero_indices[1]
                z_coords = nonzero_indices[2]
                values = cumulative[x_coords, y_coords, z_coords]
            
            # 创建3D散点图
            ax = plt.gca(projection='3d')
            scatter = ax.scatter(x_coords, y_coords, z_coords, c=values, 
                               cmap='viridis', s=20, alpha=0.6)
            ax.set_xlabel('X网格')
            ax.set_ylabel('Y网格')
            ax.set_zlabel('Z网格')
            ax.set_title('3D网格数据分布')
            plt.colorbar(scatter)
    
    # X方向积分
    plt.subplot(2, 2, 4)
    if len(cumulative.shape) == 3:
        x_integrated = np.sum(cumulative, axis=(1, 2))
        plt.plot(x_integrated)
        plt.xlabel('X网格索引')
        plt.ylabel('积分值')
        plt.title('X方向积分')
        plt.grid(True)
    
    plt.tight_layout()
    plt.show()
    
    # 显示详细统计信息
    print(f"\n详细网格统计:")
    print(f"  网格总体积: {np.prod(shape)}")
    print(f"  有数据的网格点: {np.count_nonzero(cumulative)}")
    print(f"  数据覆盖率: {np.count_nonzero(cumulative)/np.prod(shape)*100:.2f}%")
    
    if len(cumulative.shape) == 3:
        print(f"\n各方向统计:")
        print(f"  X方向最大值位置: {np.argmax(np.sum(cumulative, axis=(1,2)))}")
        print(f"  Y方向最大值位置: {np.argmax(np.sum(cumulative, axis=(0,2)))}")
        print(f"  Z方向最大值位置: {np.argmax(np.sum(cumulative, axis=(0,1)))}")