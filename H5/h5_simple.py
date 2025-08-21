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
    
    # 打印第100个batch（索引99）的非零数据分析
    print(f"\n第100个Batch（索引99）的非零数据分析:")
    print("-" * 100)
    batch_100_data = f['batch_99'][:]
    flat_data = batch_100_data.flatten()
    
    print(f"Batch 99数据形状: {batch_100_data.shape}")
    print(f"展平后数据长度: {len(flat_data)}")
    
    # 找到所有非零数据的索引和值
    nonzero_indices = np.nonzero(flat_data)[0]
    nonzero_data = flat_data[nonzero_indices]
    
    print(f"\n非零数据统计:")
    print(f"  总数据点: {len(flat_data)}")
    print(f"  非零数据点: {len(nonzero_data)}")
    
    if len(nonzero_data) > 0:
        print(f"  非零数据范围:")
        print(f"    第一个非零数据位置: 索引 {nonzero_indices[0]}")
        print(f"    最后一个非零数据位置: 索引 {nonzero_indices[-1]}")
        print(f"    非零数据跨度: {nonzero_indices[-1] - nonzero_indices[0] + 1} 个位置")
        print(f"    数据密度: {len(nonzero_data)/(nonzero_indices[-1] - nonzero_indices[0] + 1)*100:.2f}%")
        
        print(f"\n前200个非零数据:")
        print("-" * 100)
        print(f"{'序号':>4} {'索引':>8} {'数值':>15} {'3D坐标(ix,iy,iz)':>20}")
        print("-" * 100)
        
        for i in range(min(200, len(nonzero_data))):
            idx = nonzero_indices[i]
            value = nonzero_data[i]
            
            # 将1D索引转换为3D坐标 (ix, iy, iz)
            iz = idx // (shape[0] * shape[1])
            remainder = idx % (shape[0] * shape[1])
            iy = remainder // shape[0]
            ix = remainder % shape[0]
            
            print(f"{i+1:>4} {idx:>8} {value:>15.3e} ({ix:>2},{iy:>2},{iz:>2})")
        
        print("-" * 100)
        print(f"显示了前 {min(200, len(nonzero_data))} 个非零数据（共 {len(nonzero_data)} 个）")
        
        # 按Z方向分层显示2D矩阵
        print(f"\n按Z方向分层的2D矩阵显示:")
        print("=" * 120)
        
        # 重新整形数据为3D
        data_3d = batch_100_data.reshape(shape[0], shape[1], shape[2])
        
        # 找到有数据的Z层
        z_layers_with_data = []
        for z in range(shape[2]):
            layer_data = data_3d[:, :, z]
            if np.any(layer_data > 0):
                z_layers_with_data.append(z)
        
        print(f"发现 {len(z_layers_with_data)} 个Z层有数据: {z_layers_with_data}")
        
        for z in z_layers_with_data:
            layer_data = data_3d[:, :, z]
            nonzero_count = np.count_nonzero(layer_data)
            layer_sum = np.sum(layer_data)
            layer_max = np.max(layer_data)
            
            print(f"\nZ层 {z} (物理坐标: {origin[2] + z * f.attrs['pitch']:.2f} cm)")
            print(f"非零点数: {nonzero_count}, 总和: {layer_sum:.2e}, 最大值: {layer_max:.2e}")
            print("-" * 120)
            
            # 找到有数据的X和Y范围
            nonzero_x, nonzero_y = np.nonzero(layer_data)
            if len(nonzero_x) > 0:
                x_min, x_max = nonzero_x.min(), nonzero_x.max()
                y_min, y_max = nonzero_y.min(), nonzero_y.max()
                
                # 扩展显示范围（加上边界）
                x_start = max(0, x_min - 2)
                x_end = min(shape[0], x_max + 3)
                y_start = max(0, y_min - 2)
                y_end = min(shape[1], y_max + 3)
                
                print(f"显示范围: X[{x_start}:{x_end}] Y[{y_start}:{y_end}]")
                
                # 打印Y坐标标题
                print("    ", end="")
                for y in range(y_start, y_end):
                    print(f"{y:>8}", end="")
                print()
                
                # 打印矩阵数据
                for x in range(x_start, x_end):
                    print(f"{x:>3} ", end="")
                    for y in range(y_start, y_end):
                        value = layer_data[x, y]
                        if value > 0:
                            print(f"{value:>8.1f}", end="")
                        else:
                            print(f"{'·':>8}", end="")  # 用点表示零值
                    print()
                
                print()
    else:
        print("  没有非零数据")

    # 