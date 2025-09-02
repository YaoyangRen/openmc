import h5py
import numpy as np
import pandas as pd

# ====== 用户可修改参数 ======
batch_idx = 99  # 第100个batch（索引从0开始）
z_indices = list(range(20))  # 从0到19
# ==========================

with h5py.File('green_function_data.h5', 'r') as f:
    shape = f['shape'][:]
    origin = f['origin'][:]
    pitch = f.attrs['pitch']
    batch_data = f[f'batch_{batch_idx}'][:]
    data_3d = batch_data.reshape(shape[0], shape[1], shape[2])
    
    # 创建输出文件名
    excel_filename = f"batch_{batch_idx}_multi_layers.xlsx"
    print(f"将处理 {len(z_indices)} 个Z层，保存到文件: {excel_filename}")
    
    # 输出全部范围（0到70）
    x_start, x_end = 0, shape[0]
    y_start, y_end = 0, shape[1]
    
    # 计算物理坐标用于列名
    y_coords = [origin[1] + y * pitch for y in range(y_start, y_end)]
    column_names = ['X\\Y_phys'] + [f"{y_coord:.1f}" for y_coord in y_coords]
    
    # 准备汇总信息
    summary_data = []
    
    # 保存到Excel文件
    with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
        for i, z_idx in enumerate(z_indices):
            layer = data_3d[:, :, z_idx]
            z_coord = origin[2] + z_idx * pitch
            nonzero_count = np.count_nonzero(layer)
            total_sum = np.sum(layer)
            max_val = np.max(layer)
            
            print(f"处理Z层 {z_idx} (物理坐标: {z_coord:.2f} cm), 非零点数: {nonzero_count}")
            
            # 创建DataFrame用于当前Z层
            df_data = []
            for x in range(x_start, x_end):
                # 计算X的物理坐标
                x_coord = origin[0] + x * pitch
                row = [f"{x_coord:.1f}"]  # 第一列是X的物理坐标
                for y in range(y_start, y_end):
                    v = layer[x, y]
                    if v > 0:
                        row.append(v)
                    else:
                        row.append(0)  # Excel中用0表示
                df_data.append(row)
            
            # 创建当前Z层的DataFrame
            df = pd.DataFrame(df_data, columns=column_names)
            
            # 保存到对应的工作表
            sheet_name = f"Z{z_idx}_coord{z_coord:.1f}cm"
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            # 添加到汇总信息
            summary_data.append({
                'Z层索引': z_idx,
                '物理坐标(cm)': f"{z_coord:.2f}",
                '非零点数': nonzero_count,
                '总和': f"{total_sum:.2e}",
                '最大值': f"{max_val:.2e}",
                '工作表名': sheet_name
            })
        
        # 创建汇总信息表
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # 创建整体信息表
        info_data = {
            '参数': ['Batch编号', 'X网格数', 'Y网格数', 'Z网格数', 'Z层数量', '网格间距'],
            '值': [batch_idx, shape[0], shape[1], shape[2], len(z_indices), pitch]
        }
        info_df = pd.DataFrame(info_data)
        info_df.to_excel(writer, sheet_name='Info', index=False)
    
    print(f"数据已成功保存到Excel文件: {excel_filename}")
    print(f"包含 {len(z_indices)} 个Z层工作表 + Summary + Info 工作表")
