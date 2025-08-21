import h5py
import numpy as np
import matplotlib.pyplot as plt

print("正在加载数据...")

with h5py.File('green_function_data.h5', 'r') as f:
    # 读取网格参数
    origin = f['origin'][:]
    shape = f['shape'][:]
    batch_100_data = f['batch_99'][:]
    
    print(f"网格形状: {shape}")
    
    # 重新整形数据为3D
    data_3d = batch_100_data.reshape(shape[0], shape[1], shape[2])
    
    # 找到非零数据最多的Z层
    print("查找非零数据最多的Z层...")
    max_count = 0
    best_z = 0
    
    for z in range(shape[2]):
        layer_data = data_3d[:, :, z]
        nonzero_count = np.count_nonzero(layer_data)
        if nonzero_count > max_count:
            max_count = nonzero_count
            best_z = z
    
    # 获取最佳Z层的数据
    layer_data = data_3d[:, :, best_z]
    layer_sum = np.sum(layer_data)
    layer_max = np.max(layer_data)
    
    print(f"选择Z层 {best_z}:")
    print(f"  物理坐标: {origin[2] + best_z * f.attrs['pitch']:.1f} cm")
    print(f"  非零点数: {max_count}")
    print(f"  总和: {layer_sum:.1f}")
    print(f"  最大值: {layer_max:.1f}")
    
    # 找到有数据的区域
    nonzero_x, nonzero_y = np.nonzero(layer_data)
    x_min, x_max = nonzero_x.min(), nonzero_x.max()
    y_min, y_max = nonzero_y.min(), nonzero_y.max()
    
    print(f"数据分布范围: X[{x_min}:{x_max}], Y[{y_min}:{y_max}]")
    
    # 扩展显示范围
    x_start = max(0, x_min - 3)
    x_end = min(shape[0], x_max + 4)
    y_start = max(0, y_min - 3)
    y_end = min(shape[1], y_max + 4)
    
    print(f"显示范围: X[{x_start}:{x_end}], Y[{y_start}:{y_end}]")
    
    # 提取显示区域的数据
    display_data = layer_data[x_start:x_end, y_start:y_end]
    
    print("正在生成热力图...")
    
    # 创建热力图
    plt.figure(figsize=(14, 10))
    
    # 使用热力图颜色方案
    im = plt.imshow(display_data, cmap='hot', origin='lower', aspect='equal', interpolation='nearest')
    
    # 添加颜色条
    cbar = plt.colorbar(im, shrink=0.8)
    cbar.set_label('格林函数值', fontsize=14, fontfamily='SimHei')
    
    # 设置标题
    plt.title(f'格林函数热力图 - Z层 {best_z}\n'
             f'物理坐标: {origin[2] + best_z * f.attrs["pitch"]:.1f} cm, '
             f'非零点: {max_count}个, 最大值: {layer_max:.1f}', 
             fontsize=16, fontfamily='SimHei', pad=20)
    
    # 设置坐标轴标签
    plt.xlabel(f'Y网格索引', fontsize=14, fontfamily='SimHei')
    plt.ylabel(f'X网格索引', fontsize=14, fontfamily='SimHei')
    
    # 设置坐标轴刻度
    x_ticks = np.arange(0, x_end-x_start, 3)
    y_ticks = np.arange(0, y_end-y_start, 3)
    plt.xticks(y_ticks, [y_start + i for i in y_ticks])
    plt.yticks(x_ticks, [x_start + i for i in x_ticks])
    
    # 添加网格
    plt.grid(True, alpha=0.3, color='white', linewidth=0.5)
    
    # 在所有点上标注数值
    threshold = layer_max * 0.0  # 只标注大于最大值0%的点
    labeled_count = 0
    for i in range(display_data.shape[0]):
        for j in range(display_data.shape[1]):
            value = display_data[i, j]
            if value > threshold and labeled_count < 15:  # 限制标注数量
                plt.text(j, i, f'{value:.0f}', 
                       ha='center', va='center', 
                       color='cyan' if value > layer_max * 0.7 else 'yellow',
                       fontsize=10, weight='bold',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.7))
                labeled_count += 1
    
    plt.tight_layout()
    
    # 保存图片
    filename = f'green_function_heatmap_z{best_z}_enhanced.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"热力图已保存为: {filename}")
    
    # 显示统计信息
    print(f"\n数据统计:")
    print(f"  显示区域大小: {display_data.shape}")
    print(f"  非零值范围: {np.min(display_data[display_data>0]):.2f} - {np.max(display_data):.2f}")
    print(f"  平均值(非零): {np.mean(display_data[display_data>0]):.2f}")
    
    plt.show()

print("热力图生成完成！")
