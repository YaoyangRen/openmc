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
    
    # 读取第100个batch数据
    batch_100_data = f['batch_99'][:]
    
    print(f"\n第100个Batch按Z方向分层的2D矩阵显示:")
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
    
    for z in z_layers_with_data[:5]:  # 只显示前5层
        layer_data = data_3d[:, :, z]
        nonzero_count = np.count_nonzero(layer_data)
        layer_sum = np.sum(layer_data)
        layer_max = np.max(layer_data)
        
        print(f"\nZ层 {z} (物理坐标: {origin[2] + z * f.attrs['pitch']:.2f} cm)")
        print(f"非零点数: {nonzero_count}, 总和: {layer_sum:.2e}, 最大值: {layer_max:.2e}")
        print("-" * 80)
        
        # 找到有数据的X和Y范围
        nonzero_x, nonzero_y = np.nonzero(layer_data)
        if len(nonzero_x) > 0:
            x_min, x_max = nonzero_x.min(), nonzero_x.max()
            y_min, y_max = nonzero_y.min(), nonzero_y.max()
            
            # 限制显示范围
            x_start = max(0, x_min - 1)
            x_end = min(shape[0], x_max + 2)
            y_start = max(0, y_min - 1)
            y_end = min(shape[1], y_max + 2)
            
            # 进一步限制显示大小
            if (x_end - x_start) > 20:
                x_end = x_start + 20
            if (y_end - y_start) > 15:
                y_end = y_start + 15
                
            print(f"显示范围: X[{x_start}:{x_end}] Y[{y_start}:{y_end}]")
            
            # 打印Y坐标标题
            print("    ", end="")
            for y in range(y_start, y_end):
                print(f"{y:>6}", end="")
            print()
            
            # 打印矩阵数据
            for x in range(x_start, x_end):
                print(f"{x:>3} ", end="")
                for y in range(y_start, y_end):
                    value = layer_data[x, y]
                    if value > 0:
                        if value >= 1000:
                            print(f"{value:>6.0f}", end="")
                        elif value >= 1:
                            print(f"{value:>6.1f}", end="")
                        else:
                            print(f"{value:>6.2f}", end="")
                    else:
                        print(f"{'·':>6}", end="")  # 用点表示零值
                print()
            
            print()
    
    if len(z_layers_with_data) > 5:
        print(f"... 还有 {len(z_layers_with_data)-5} 个Z层有数据")
