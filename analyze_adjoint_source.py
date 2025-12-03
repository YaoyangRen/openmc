"""
伴随源迭代计算示例
===================

演示如何使用裂变矩阵计算伴随源分布

理论基础：
    伴随裂变源满足：I* = (1/k) F^T I*
    
    其中：
    - I* 是伴随源分布向量
    - F^T 是裂变矩阵的转置
    - k 是有效增殖因子

初始值设置：
    1. 均匀分布 (uniform): I*_i = 1
    2. 正向源分布 (forward): I*_i = S_i
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix

def read_fission_matrix(filename='fission_matrix.h5'):
    """读取裂变矩阵和伴随源数据"""
    
    with h5py.File(filename, 'r') as f:
        # 读取网格信息
        origin = f['origin'][:]
        shape = f['shape'][:]
        pitch = f.attrs['pitch']
        n_cells = f.attrs['n_cells']
        nnz = f.attrs['nnz']
        
        print("=" * 70)
        print("FISSION MATRIX AND ADJOINT SOURCE DATA")
        print("=" * 70)
        print(f"\nGrid Information:")
        print(f"  Shape: {shape}")
        print(f"  Origin: {origin}")
        print(f"  Pitch: {pitch} cm")
        print(f"  Total cells: {n_cells}")
        print(f"  Non-zero elements: {nnz}")
        
        # 读取稀疏矩阵 (COO格式)
        row_indices = f['row_indices'][:]
        col_indices = f['col_indices'][:]
        data_normalized = f['data_normalized'][:]
        
        # 构建稀疏矩阵
        F = coo_matrix(
            (data_normalized, (row_indices, col_indices)),
            shape=(n_cells, n_cells)
        )
        
        # 读取源计数
        source_counts = f['source_counts'][:]
        
        # 读取伴随源（如果存在）
        adjoint_source = None
        keff_reference = None
        adjoint_iterations = None
        adjoint_history = None
        adjoint_history_label = None
        batch_adjoint_enabled = False
        adjoint_iter_per_batch = 0
        
        if 'adjoint_source' in f:
            adjoint_source = f['adjoint_source'][:]
            keff_reference = f.attrs.get('reference_keff')
            if keff_reference is None:
                keff_reference = f.attrs.get('keff_reference')
            if keff_reference is None:
                keff_reference = f.attrs.get('k_adjoint')  # backward compat
            adjoint_iterations = f.attrs.get('adjoint_iterations')
            
            print(f"\nAdjoint Source Information:")
            if keff_reference is not None:
                print(f"  keff_reference = {keff_reference:.8f}")
            if adjoint_iterations is not None:
                print(f"  Iterations: {adjoint_iterations}")
            else:
                print("  Iterations: N/A")
            print(f"  Adjoint source available: Yes")
            
            # 读取batch级迭代历史（如果存在）
            if 'adjoint_residual_history' in f:
                adjoint_history = f['adjoint_residual_history'][:] # type: ignore
                adjoint_history_label = 'residual'
            elif 'k_adjoint_history' in f:
                adjoint_history = f['k_adjoint_history'][:]  # type: ignore[index]
                adjoint_history_label = 'legacy_k'
                batch_adjoint_enabled = f.attrs.get('batch_adjoint_enabled', False)
                adjoint_iter_per_batch = f.attrs.get('adjoint_iter_per_batch', 0)
                
                print(f"\nBatch-level Adjoint Iteration:")
                print(f"  Enabled: {batch_adjoint_enabled}")
                print(f"  Iterations per batch: {adjoint_iter_per_batch}")
                print(f"  Number of batches: {len(adjoint_history)}") # type: ignore
                if adjoint_history_label == 'residual':
                    print(f"  Residual evolution (max |ΔI*| per batch):")
                    print(f"    Initial = {adjoint_history[0]:.3e}")
                    print(f"    Final   = {adjoint_history[-1]:.3e}")
                elif adjoint_history_label == 'legacy_k':
                    print(f"  Legacy k_adjoint evolution:")
                    print(f"    Initial = {adjoint_history[0]:.8f}")
                    print(f"    Final   = {adjoint_history[-1]:.8f}")
                    print(f"    Change  = {abs(adjoint_history[-1] - adjoint_history[0]):.3e}")
        else:
            print(f"\nAdjoint source not computed")
        
        print("=" * 70)
        
        return {
            'F': F,
            'shape': shape,
            'origin': origin,
            'pitch': pitch,
            'source_counts': source_counts,
            'adjoint_source': adjoint_source,
            'keff_reference': keff_reference,
            'adjoint_iterations': adjoint_iterations,
            'adjoint_history': adjoint_history,
            'adjoint_history_label': adjoint_history_label,
            'batch_adjoint_enabled': batch_adjoint_enabled,
            'adjoint_iter_per_batch': adjoint_iter_per_batch
        }

def compute_adjoint_manually(
    F, keff, initial_guess='uniform', max_iter=1000, tol=1e-6):
    """
    手动计算伴随源（用于验证C++实现）
    
    Parameters:
    -----------
    F : scipy.sparse matrix
        裂变矩阵 (已归一化)
    initial_guess : str
        'uniform' 或 'forward'
    max_iter : int
        最大迭代次数
    tol : float
        收敛容差
    """
    
    n = F.shape[0]
    F_csr = F.tocsr()  # 转换为CSR格式以加速矩阵-向量乘法
    
    # 初始化
    if initial_guess == 'uniform':
        I_star = np.ones(n)
    elif initial_guess == 'forward':
        # 使用F的列和作为正向源
        I_star = np.array(F.sum(axis=0)).flatten()
        if I_star.sum() == 0:
            I_star = np.ones(n)
    else:
        I_star = np.ones(n)
    
    # 归一化
    I_star = I_star / I_star.sum()
    
    print(f"\nManual Adjoint Computation")
    print(f"  Initial guess: {initial_guess}")
    print(f"  Max iterations: {max_iter}")
    print(f"  Tolerance: {tol:.2e}")
    print(f"  Reference keff: {keff:.8f}")
    
    max_delta = np.inf
    for iteration in range(max_iter):
        # I_new = F^T × I*
        I_new = F_csr.T.dot(I_star)
        
        I_new /= keff
        sum_new = I_new.sum()
        
        if sum_new == 0:
            print(f"Error: adjoint collapsed to zero at iteration {iteration}")
            break
        
        # 归一化
        new_vector = I_new / sum_new
        max_delta = np.max(np.abs(new_vector - I_star))
        I_star = new_vector
        
        # 检查收敛
        if (iteration + 1) % 50 == 0 or iteration == 0:
            print(f"  Iteration {iteration+1:4d}: max |ΔI*| = {max_delta:.2e}")
        
        if max_delta < tol:
            print(f"\nConverged at iteration {iteration+1}")
            print(f"  Final max |ΔI*| = {max_delta:.2e}")
            return I_star, max_delta, iteration + 1
    
    print(f"\nWarning: Maximum iterations reached")
    print(f"  Final max |ΔI*| = {max_delta:.2e}")
    return I_star, max_delta, max_iter

def visualize_adjoint_source(data, slice_axis='z', slice_index=None):
    """
    可视化伴随源分布
    
    Parameters:
    -----------
    data : dict
        包含伴随源和网格信息的字典
    slice_axis : str
        切片轴 ('x', 'y', 'z')
    slice_index : int
        切片索引（默认为中间位置）
    """
    
    if data['adjoint_source'] is None:
        print("Error: Adjoint source not available")
        return
    
    shape = data['shape']
    adjoint_source = data['adjoint_source']
    
    # 重塑为3D数组
    I_3d = adjoint_source.reshape(shape)
    
    # 确定切片索引
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_idx = axis_map[slice_axis]
    
    if slice_index is None:
        slice_index = shape[axis_idx] // 2
    
    # 提取切片
    if slice_axis == 'x':
        I_slice = I_3d[slice_index, :, :]
        xlabel, ylabel = 'Y', 'Z'
    elif slice_axis == 'y':
        I_slice = I_3d[:, slice_index, :]
        xlabel, ylabel = 'X', 'Z'
    else:  # z
        I_slice = I_3d[:, :, slice_index]
        xlabel, ylabel = 'X', 'Y'
    
    # 绘图
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 2D热图
    im1 = axes[0].imshow(I_slice.T, origin='lower', cmap='hot', 
                         interpolation='nearest')
    axes[0].set_xlabel(f'{xlabel} index')
    axes[0].set_ylabel(f'{ylabel} index')
    axes[0].set_title(f'Adjoint Source ({slice_axis.upper()}-slice at index {slice_index})')
    plt.colorbar(im1, ax=axes[0], label='Adjoint Importance')
    
    # 对数刻度热图（突出小值）
    I_slice_log = np.log10(I_slice + 1e-20)  # 避免log(0)
    im2 = axes[1].imshow(I_slice_log.T, origin='lower', cmap='viridis',
                         interpolation='nearest')
    axes[1].set_xlabel(f'{xlabel} index')
    axes[1].set_ylabel(f'{ylabel} index')
    axes[1].set_title(f'Log10(Adjoint Source)')
    plt.colorbar(im2, ax=axes[1], label='Log10(Importance)')
    
    keff_value = data.get('keff_reference')
    iteration_info = data.get('adjoint_iterations')
    title_parts = []
    if keff_value is not None:
        title_parts.append(f'keff_reference = {keff_value:.6f}')
    if iteration_info is not None:
        title_parts.append(f'Iterations = {iteration_info}')
    title_text = ', '.join(title_parts) if title_parts else 'Adjoint Source'
    plt.suptitle(title_text, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('adjoint_source_distribution.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved as 'adjoint_source_distribution.png'")
    plt.show()

def compare_forward_adjoint(data):
    """比较正向源和伴随源分布"""
    
    if data['adjoint_source'] is None:
        print("Error: Adjoint source not available")
        return
    
    shape = data['shape']
    source_counts = data['source_counts']
    adjoint_source = data['adjoint_source']
    
    # 归一化正向源
    forward_source = source_counts / source_counts.sum() if source_counts.sum() > 0 else source_counts
    
    # 重塑为3D
    F_3d = forward_source.reshape(shape)
    I_3d = adjoint_source.reshape(shape)
    
    # 提取中间切片
    mid_z = shape[2] // 2
    F_slice = F_3d[:, :, mid_z]
    I_slice = I_3d[:, :, mid_z]
    
    # 绘图对比
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # 正向源
    im1 = axes[0].imshow(F_slice.T, origin='lower', cmap='Reds', 
                         interpolation='nearest')
    axes[0].set_title('Forward Source Distribution')
    axes[0].set_xlabel('X index')
    axes[0].set_ylabel('Y index')
    plt.colorbar(im1, ax=axes[0], label='Forward Importance')
    
    # 伴随源
    im2 = axes[1].imshow(I_slice.T, origin='lower', cmap='Blues',
                         interpolation='nearest')
    axes[1].set_title('Adjoint Source Distribution')
    axes[1].set_xlabel('X index')
    axes[1].set_ylabel('Y index')
    plt.colorbar(im2, ax=axes[1], label='Adjoint Importance')
    
    # 比值 (I*/S)
    ratio = np.divide(I_slice, F_slice, 
                      out=np.zeros_like(I_slice), 
                      where=F_slice!=0)
    im3 = axes[2].imshow(ratio.T, origin='lower', cmap='RdYlGn',
                         interpolation='nearest')
    axes[2].set_title('Adjoint / Forward Ratio')
    axes[2].set_xlabel('X index')
    axes[2].set_ylabel('Y index')
    plt.colorbar(im3, ax=axes[2], label='I* / S')
    
    plt.suptitle(f'Forward vs Adjoint Source (Z-slice at {mid_z})',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('forward_adjoint_comparison.png', dpi=150, bbox_inches='tight')
    print("Comparison plot saved as 'forward_adjoint_comparison.png'")
    plt.show()

def analyze_convergence_history(data):
    """分析伴随源迭代收敛历史（需要修改C++代码保存历史）"""
    
    # 这里可以手动重新计算以获得收敛历史
    F = data['F']
    keff = data.get('keff_reference') or 1.0
    
    print("\nRecomputing adjoint for convergence analysis...")
    
    # 使用均匀初始化
    I_uniform, res_uniform, iter_uniform = compute_adjoint_manually(
        F, keff, initial_guess='uniform', max_iter=1000, tol=1e-8
    )
    
    print("\nUsing forward source initialization...")
    I_forward, res_forward, iter_forward = compute_adjoint_manually(
        F, keff, initial_guess='forward', max_iter=1000, tol=1e-8
    )
    
    print(f"\nComparison:")
    print(f"  Uniform init: {iter_uniform} iterations, max |ΔI*| = {res_uniform:.2e}")
    print(f"  Forward init: {iter_forward} iterations, max |ΔI*| = {res_forward:.2e}")
    print(f"  L2 norm of I* difference: {np.linalg.norm(I_uniform - I_forward):.2e}")

def plot_batch_adjoint_convergence(data):
    """
    绘制batch级伴随源迭代的收敛历史
    
    Parameters:
    -----------
    data : dict
        从read_fission_matrix()返回的数据字典
    """
    
    history = data.get('adjoint_history')
    history_label = data.get('adjoint_history_label')
    
    if history is None:
        print("\nNo batch-level adjoint iteration history available.")
        print("To enable, call in C++:")
        print("  fission_matrix->enable_batch_adjoint_iteration(true, 10, 1e-6);")
        return
    
    print("\nPlotting batch-level adjoint convergence...")
    
    n_batches = len(history)
    batches = np.arange(1, n_batches + 1)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    if history_label == 'residual':
        ax1.semilogy(batches, history, 'b-o', linewidth=2, markersize=4)
        ax1.set_ylabel('max |ΔI*|', fontsize=12)
        ax1.set_title('Batch-wise Adjoint Residual', fontsize=14)
        ax1.grid(True, alpha=0.3, which='both')

        ax2.plot(batches, history, 'g-o', linewidth=2, markersize=4)
        ax2.set_xlabel('Batch Number', fontsize=12)
        ax2.set_ylabel('max |ΔI*|', fontsize=12)
        ax2.set_title('Residual (Linear Scale)', fontsize=14)
        ax2.grid(True, alpha=0.3)
    else:
        ax1.plot(batches, history, 'b-o', linewidth=2, markersize=4)
        ax1.set_xlabel('Batch Number', fontsize=12)
        ax1.set_ylabel('$k_{adjoint}$', fontsize=12)
        ax1.set_title('Adjoint Eigenvalue (Legacy)', fontsize=14)
        ax1.grid(True, alpha=0.3)

        final_k = history[-1]
        ax1.axhline(y=final_k, color='r', linestyle='--', alpha=0.5,
                    label=f'Final: {final_k:.8f}')
        ax1.legend()

        if n_batches > 1:
            relative_error = np.abs((history - final_k) / final_k)
            ax2.semilogy(
                batches, relative_error, 'g-o', linewidth=2, markersize=4)
            ax2.set_xlabel('Batch Number', fontsize=12)
            ax2.set_ylabel('Relative Error $|k - k_{final}| / k_{final}$', fontsize=12)
            ax2.set_title('Convergence Rate', fontsize=14)
            ax2.grid(True, alpha=0.3, which='both')

            tolerance = 1e-6
            ax2.axhline(y=tolerance, color='r', linestyle='--', alpha=0.5,
                        label=f'Tolerance: {tolerance:.1e}')
            ax2.legend()
    
    plt.tight_layout()
    plt.savefig('batch_adjoint_convergence.png', dpi=300, bbox_inches='tight')
    print(f"  Saved: batch_adjoint_convergence.png")
    plt.close()
    
    # 打印统计信息
    print(f"\nBatch-level Convergence Statistics:")
    print(f"  Number of batches: {n_batches}")
    if history_label == 'residual':
        print(f"  Initial residual: {history[0]:.3e}")
        print(f"  Final residual:   {history[-1]:.3e}")
    else:
        print(f"  Initial k_adjoint: {history[0]:.8f}")
        print(f"  Final k_adjoint:   {history[-1]:.8f}")
        print(f"  Total change:      {abs(history[-1] - history[0]):.3e}")

        if n_batches > 10:
            last_10 = history[-10:]
            std_last_10 = np.std(last_10)
            mean_last_10 = np.mean(last_10)
            rel_std = std_last_10 / mean_last_10
            print(f"\nLast 10 batches stability:")
            print(f"  Mean:     {mean_last_10:.8f}")
            print(f"  Std dev:  {std_last_10:.3e}")
            print(f"  Rel std:  {rel_std:.3e}")

# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    import sys
    
    # 读取数据
    filename = 'fission_matrix.h5'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    
    print(f"\nReading fission matrix from: {filename}")
    data = read_fission_matrix(filename)
    
    if data['adjoint_source'] is not None:
        print("\n" + "=" * 70)
        print("VISUALIZATION")
        print("=" * 70)
        
        # 可视化伴随源分布
        visualize_adjoint_source(data, slice_axis='z')
        
        # 比较正向和伴随源
        compare_forward_adjoint(data)
        
        # 绘制batch级迭代收敛历史（如果有）
        if data.get('adjoint_history') is not None:
            plot_batch_adjoint_convergence(data)
        
        # 分析收敛性（可选）
        print("\n" + "=" * 70)
        print("CONVERGENCE ANALYSIS")
        print("=" * 70)
        analyze_convergence_history(data)
    else:
        print("\nAdjoint source has not been computed.")
        print("Please run compute_adjoint_source() in C++ before calling finalize().")
        print("\nExample C++ usage:")
        print("  fission_matrix->compute_adjoint_source(\"uniform\", 1000, 1e-6);")
        print("  fission_matrix->finalize(\"fission_matrix.h5\");")
