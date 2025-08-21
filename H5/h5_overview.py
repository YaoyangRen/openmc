import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File('green_function_data.h5', 'r') as f:
    # 读取网格参数
    origin = f['origin'][:]
    shape = f['shape'][:]
    batch_100_data = f['batch_99'][:]
    
    print(f"第100个Batch的网格数据概览 (网格形状: {shape})")
    print("=" * 100)
    
    # 重新整形数据为3D
    data_3d = batch_100_data.reshape(shape[0], shape[1], shape[2])
    
    # 找到有数据的Z层并按数据量排序
    z_layer_stats = []
    for z in range(shape[2]):
        layer_data = data_3d[:, :, z]
        nonzero_count = np.count_nonzero(layer_data)
        if nonzero_count > 0:
            layer_sum = np.sum(layer_data)
            layer_max = np.max(layer_data)
            z_layer_stats.append((z, nonzero_count, layer_sum, layer_max))
    
    # 按数据点数量排序，显示数据最多的层
    z_layer_stats.sort(key=lambda x: x[1], reverse=True)
    
    print(f"发现 {len(z_layer_stats)} 个Z层有数据")
    print(f"数据最多的前10个Z层:")
    for i, (z, count, sum_val, max_val) in enumerate(z_layer_stats[:10]):
        print(f"  Z{z:>2} (物理:{origin[2] + z * f.attrs['pitch']:>6.1f}cm): {count:>3}个点, 总和:{sum_val:>8.1f}, 最大:{max_val:>8.1f}")
    
    # 显示数据最集中的几层的详细矩阵
    print(f"\n详细矩阵显示 (前3个数据最多的Z层):")
    
    for i, (z, count, sum_val, max_val) in enumerate(z_layer_stats[:3]):
        layer_data = data_3d[:, :, z]
        
        print(f"\nZ层 {z} (物理坐标: {origin[2] + z * f.attrs['pitch']:.1f} cm)")
        print(f"非零点数: {count}, 总和: {sum_val:.1f}, 最大值: {max_val:.1f}")
        print("-" * 80)
        
        # 找到有数据的区域
        nonzero_x, nonzero_y = np.nonzero(layer_data)
        x_min, x_max = nonzero_x.min(), nonzero_x.max()
        y_min, y_max = nonzero_y.min(), nonzero_y.max()
        
        # 设置显示范围
        x_start = max(0, x_min - 1)
        x_end = min(shape[0], x_max + 2)
        y_start = max(0, y_min - 1)  
        y_end = min(shape[1], y_max + 2)
        
        # 限制显示大小
        if (x_end - x_start) > 15:
            x_center = (x_min + x_max) // 2
            x_start = max(0, x_center - 7)
            x_end = min(shape[0], x_center + 8)
        if (y_end - y_start) > 20:
            y_center = (y_min + y_max) // 2
            y_start = max(0, y_center - 10)
            y_end = min(shape[1], y_center + 11)
            
        print(f"显示区域中心: X[{x_start}:{x_end}] Y[{y_start}:{y_end}]")
        
        # 打印Y坐标
        print("   ", end="")
        for y in range(y_start, y_end):
            print(f"{y:>5}", end="")
        print()
        
        # 打印矩阵
        for x in range(x_start, x_end):
            print(f"{x:>2} ", end="")
            for y in range(y_start, y_end):
                value = layer_data[x, y]
                if value > 0:
                    if value >= 1000:
                        print(f"{value:>5.0f}", end="")
                    elif value >= 10:
                        print(f"{value:>5.1f}", end="")
                    else:
                        print(f"{value:>5.2f}", end="")
                else:
                    print(f"{'·':>5}", end="")
            print()
        print()
    
    # 显示Z方向分布概览
    print(f"\nZ方向数据分布概览:")
    print("Z层范围     非零点数  总数据量")
    print("-" * 40)
    
    # 按Z层分组显示
    for start_z in range(0, len(z_layer_stats), 10):
        end_z = min(start_z + 10, len(z_layer_stats))
        group_layers = z_layer_stats[start_z:end_z]
        
        if group_layers:
            min_z = min(layer[0] for layer in group_layers)
            max_z = max(layer[0] for layer in group_layers)
            total_points = sum(layer[1] for layer in group_layers)
            total_sum = sum(layer[2] for layer in group_layers)
            
            print(f"Z{min_z:>2}-{max_z:<2}      {total_points:>8}  {total_sum:>8.1f}")
    
    # 为非零数据最多的Z层创建热力图
    if z_layer_stats:
        best_z, best_count, best_sum, best_max = z_layer_stats[0]
        layer_data = data_3d[:, :, best_z]
        
        print(f"\n为Z层 {best_z} 创建热力图...")
        
        # 创建热力图
        plt.figure(figsize=(12, 10))
        
        # 找到有数据的区域来设置合适的显示范围
        nonzero_x, nonzero_y = np.nonzero(layer_data)
        if len(nonzero_x) > 0:
            x_min, x_max = nonzero_x.min(), nonzero_x.max()
            y_min, y_max = nonzero_y.min(), nonzero_y.max()
            
            # 扩展显示范围
            x_start = max(0, x_min - 5)
            x_end = min(shape[0], x_max + 6)
            y_start = max(0, y_min - 5)
            y_end = min(shape[1], y_max + 6)
            
            # 提取显示区域的数据
            display_data = layer_data[x_start:x_end, y_start:y_end]
            
            # 创建热力图
            im = plt.imshow(display_data, cmap='hot', origin='lower', aspect='equal')
            
            # 添加颜色条
            cbar = plt.colorbar(im)
            cbar.set_label('绿函数值', fontsize=12)
            
            # 设置标题和标签
            plt.title(f'绿函数热力图 - Z层 {best_z} (物理坐标: {origin[2] + best_z * f.attrs["pitch"]:.1f} cm)\n'
                     f'非零点数: {best_count}, 总和: {best_sum:.1f}, 最大值: {best_max:.1f}', 
                     fontsize=14, pad=20)
            
            # 设置坐标轴标签（对应实际的网格索引）
            plt.xlabel(f'Y网格索引 (显示范围: {y_start}-{y_end-1})', fontsize=12)
            plt.ylabel(f'X网格索引 (显示范围: {x_start}-{x_end-1})', fontsize=12)
            
            # 设置坐标轴刻度
            x_ticks = np.arange(0, x_end-x_start, 2)
            y_ticks = np.arange(0, y_end-y_start, 2)
            plt.xticks(y_ticks, [y_start + i for i in y_ticks])
            plt.yticks(x_ticks, [x_start + i for i in x_ticks])
            
            # 在非零点上标注数值（只标注较大的值以避免过于拥挤）
            threshold = best_max * 0.3  # 只标注大于最大值30%的点
            for i in range(display_data.shape[0]):
                for j in range(display_data.shape[1]):
                    value = display_data[i, j]
                    if value > threshold:
                        plt.text(j, i, f'{value:.0f}', 
                               ha='center', va='center', 
                               color='white' if value > best_max * 0.7 else 'black',
                               fontsize=8, weight='bold')
            
            plt.tight_layout()
            
            # 保存图片
            plt.savefig(f'green_function_heatmap_z{best_z}.png', dpi=300, bbox_inches='tight')
            print(f"热力图已保存为: green_function_heatmap_z{best_z}.png")
            
            plt.show()
        
        else:
            print("没有找到非零数据！")
