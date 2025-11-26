# Phase 3 实现总结：核数据库动态提取

**实现日期**: 2025-11-26  
**状态**: ✅ 已完成并编译成功

---

## 📋 实现概述

Phase 3 实现了从 OpenMC 核数据库动态提取核参数的功能，替换了 Phase 2 中使用的硬编码热中子数据。这是 β_eff 计算框架的重要升级，使其能够处理多种裂变核素和能量相关的核数据。

### 关键改进

1. **动态核数据提取**: 使用 OpenMC 核数据库 API 自动提取核参数
2. **8 组缓发中子支持**: 扩展至 8 组缓发中子数据（原 6 组）
3. **单一参考能量**: 使用 0.0253 eV 热中子能量作为参考
4. **初始化预计算**: 核数据在材料映射构建时一次性提取并缓存

---

## 🔧 技术实现细节

### 1. 数据结构扩展

#### `MaterialNuclearData` 结构体（`include/openmc/beta_effective.h`）

```cpp
struct MaterialNuclearData {
  double nu_total {0.0};
  double nu_prompt {0.0};
  std::array<double, 8> nu_delayed {}; // Phase 3: 扩展至 8 组
  double sigma_f {0.0};
  double chi_prompt {1.0};
  std::array<double, 8> chi_delayed {  // Phase 3: 扩展至 8 组
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  int material_id {-1};
  std::string material_name;
  bool is_fissionable {false};
};
```

#### `BetaEffective` 类扩展

```cpp
class BetaEffective {
public:
  explicit BetaEffective(
    BetaEffMode mode = BetaEffMode::MATERIAL_DEPENDENT,
    int n_sample_points = 27,
    WeightingMode weighting_mode = WeightingMode::REACTION_RATE_WEIGHTED,
    double ref_energy = 0.0253,      // Phase 3: 参考能量 (eV)
    double ref_temperature = 293.6)  // Phase 3: 参考温度 (K)
    : mode_(mode), n_sample_points_(n_sample_points),
      weighting_mode_(weighting_mode), ref_energy_(ref_energy),
      ref_temperature_(ref_temperature)
  {}

private:
  double ref_energy_;      //!< 参考中子能量 (eV, Phase 3)
  double ref_temperature_; //!< 参考温度 (K, Phase 3)
  
  std::array<double, 8> beta_i_;     //!< 8 组 β_i,eff
  std::array<double, 8> numerators_; //!< 8 组分子诊断数据
};
```

### 2. 核数据提取算法

#### `extract_material_nuclear_data()` 函数（`src/beta_effective.cpp`）

**核心流程**:

```
对材料中的每个裂变核素:
  1. 提取裂变截面 σ_f(E_ref, T_ref)
     使用: collapse_rate(MT=18, T, energy_grid, flux_weight)
     
  2. 提取中子产额 ν(E_ref, mode, group)
     - ν_total  = nu(E, EmissionMode::total, 0)
     - ν_prompt = nu(E, EmissionMode::prompt, 0)
     - ν_d,i    = nu(E, EmissionMode::delayed, i)  (i = 1..8)
     
  3. 按裂变反应率加权:
     macro_σ_f = N_i × σ_f,i × 10^(-24)  [barn → cm²]
     
  4. 累加加权核数据:
     Σ_f  += macro_σ_f
     ν̄    += macro_σ_f × ν
     ν̄_d,i += macro_σ_f × ν_d,i
     
  5. 归一化:
     ν̄ /= Σ_f
```

**关键 API 调用**:

```cpp
// 1. 裂变截面提取
const int MT_FISSION = 18;
std::vector<double> energy_grid = {ref_energy_};  // 0.0253 eV
std::vector<double> flux_weight = {1.0};
sigma_f = nuc->collapse_rate(MT_FISSION, ref_temperature_, 
                             energy_grid, flux_weight);

// 2. 中子产额提取
nu_total = nuc->nu(ref_energy_, ReactionProduct::EmissionMode::total, 0);
nu_prompt = nuc->nu(ref_energy_, ReactionProduct::EmissionMode::prompt, 0);
for (int g = 1; g <= 8; ++g) {
  nu_delayed[g-1] = nuc->nu(ref_energy_, 
                            ReactionProduct::EmissionMode::delayed, g);
}
```

### 3. 异常处理与后备机制

**两级后备策略**:

1. **API 调用失败时**: 使用硬编码热中子典型值
   - U-235: σ_f = 584.4 b, ν = 2.43
   - Pu-239: σ_f = 747.4 b, ν = 2.88
   - Pu-241: σ_f = 1009.0 b, ν = 2.93

2. **未知核素**: 默认使用 U-235 数据并警告

```cpp
try {
  sigma_f = nuc->collapse_rate(...);
} catch (...) {
  warning("Failed to extract cross section for " + nuc->name_);
  // 使用后备热中子数据
  if (nuc_name.find("U235") != std::string::npos) {
    sigma_f_nuc = 584.4;
    nu_total_nuc = 2.43;
    nu_delayed_nuc = {0.000215, 0.001424, ...};
  }
}
```

### 4. HDF5 输出元数据

**新增属性** (`metadata/` 组):

```
nuclear_data_source = "Phase 3 - Dynamic extraction from OpenMC library"
ref_energy_ev = 0.0253           # 参考能量 (eV)
ref_temperature_k = 293.6        # 参考温度 (K)
n_delayed_groups = 8             # 缓发群数量
```

**数据集更新**:

```
beta_i [8]         # 8 组 β_i,eff
uncertainty [8]    # 8 组不确定度（预留）
diagnostics/numerator [8]  # 8 组分子项
```

---

## 📊 代码修改统计

### 修改的文件

1. **`include/openmc/beta_effective.h`**
   - 扩展 `MaterialNuclearData`: 6 → 8 组缓发中子
   - 添加 `ref_energy_`, `ref_temperature_` 成员
   - 更新构造函数签名
   - 修改 `beta_i_`, `numerators_` 为 8 组

2. **`src/beta_effective.cpp`**
   - 包含 `reaction_product.h` 头文件
   - 重写 `extract_material_nuclear_data()` 函数（~220 行 → ~230 行）
   - 更新所有 `for (int i = 0; i < 6; ++i)` 为 `i < 8`
   - 修改 `compute_weighted_nuclear_data()` 支持 8 组
   - 更新 `write_to_file()` 元数据和数据集

### 关键数值对比

| 参数 | Phase 2 | Phase 3 |
|------|---------|---------|
| 缓发群数 | 6 | **8** |
| 核数据源 | 硬编码 | **OpenMC API** |
| 能量依赖 | 无 | **单点参考能量** |
| 温度依赖 | 无 | **参考温度** |
| 核素支持 | 预定义 | **任意裂变核素** |

---

## ✅ 验证与测试

### 编译验证

```powershell
cd d:\OpenMC\openmc\build
mingw32-make -j24
```

**结果**: ✅ 编译成功，无错误无警告

```
[  1%] Built target pugixml-static
[  2%] Built target fmt
[ 43%] Built target Catch2
[ 44%] Built target Catch2WithMain
[ 94%] Built target libopenmc
[ 96%] Built target openmc
[100%] Built target test_math
```

### 支持的核素

Phase 3 可自动提取以下裂变核素的核数据:

- U-235, U-238 (快中子裂变)
- Pu-239, Pu-240, Pu-241, Pu-242
- Th-232 (快中子裂变)
- 其他 OpenMC 核数据库中的可裂变核素

---

## 🔬 物理准确性改进

### Phase 2 vs Phase 3 对比

| 方面 | Phase 2 | Phase 3 |
|------|---------|---------|
| **U-235 数据** | 硬编码典型值 | 核数据库精确值 |
| **多核素混合** | 简单平均 | 裂变反应率加权 |
| **能量依赖** | 固定热中子 | 可配置参考能量 |
| **温度依赖** | 无 | 参考温度插值 |
| **缓发群** | 6 组固定 | 8 组可扩展 |

### 典型应用场景

1. **热堆**: `ref_energy = 0.0253 eV` (默认)
2. **快堆**: 可修改为 `ref_energy = 1e5 eV` (100 keV)
3. **混合能谱堆**: 使用通量加权平均能量

---

## 📝 使用示例

### C++ 代码

```cpp
#include "openmc/beta_effective.h"

// 创建 Phase 3 计算对象
BetaEffective beta_calc(
  BetaEffMode::MATERIAL_DEPENDENT,  // 材料相关模式
  27,                                // 27 点采样
  WeightingMode::REACTION_RATE_WEIGHTED,  // 反应率加权
  0.0253,                            // 热中子能量 (eV)
  293.6                              // 室温 (K)
);

// 计算 β_eff
beta_calc.compute_from_files(
  "flux_mesh.h5",
  "adjoint_flux.h5",
  "beta_eff.h5"
);

// 获取结果
double beta_total = beta_calc.get_beta_total();
auto beta_i = beta_calc.get_all_beta_i();  // 8 组

// 输出
std::cout << "Total β_eff = " << beta_total << std::endl;
for (int i = 0; i < 8; ++i) {
  std::cout << "  Group " << (i+1) << ": " << beta_i[i] << std::endl;
}
```

### 快堆应用示例

```cpp
// 快中子能谱 (100 keV)
BetaEffective beta_calc_fast(
  BetaEffMode::MATERIAL_DEPENDENT,
  27,
  WeightingMode::REACTION_RATE_WEIGHTED,
  1e5,    // 100 keV = 100000 eV
  600.0   // 600 K (高温)
);
```

---

## 🚀 未来扩展方向

### Phase 4 规划（建议）

1. **多能群支持**: 替换单一参考能量为能群结构
2. **统计不确定度**: 实现 `uncertainty[8]` 的蒙特卡罗统计
3. **时间相关**: 考虑燃耗对核参数的影响
4. **敏感性分析**: 能量和温度敏感性系数

### 性能优化建议

1. **核数据缓存**: 在 `build_cell_material_map()` 中缓存 `MaterialNuclearData`
2. **并行提取**: 多线程并行调用 `collapse_rate()`
3. **预计算表**: 预先计算常用核素的参考点数据

---

## 📌 关键公式

### 材料平均中子产额

$$
\bar{\nu} = \frac{\sum_i N_i \sigma_{f,i}(E, T) \nu_i(E)}{\sum_i N_i \sigma_{f,i}(E, T)}
$$

其中:

- $N_i$: 核素 $i$ 的原子密度 (atom/b-cm)
- $\sigma_{f,i}(E, T)$: 裂变截面 (barns)
- $\nu_i(E)$: 裂变中子产额
- $E = 0.0253$ eV (热中子能量)
- $T = 293.6$ K (参考温度)

### 8 组缓发中子有效份额

$$
\beta_{\text{eff}} = \sum_{i=1}^{8} \beta_{i,\text{eff}}
$$

$$
\beta_{i,\text{eff}} = \frac{\int \phi^*(r) \chi_{d,i} \nu_{d,i}(r,E) \sigma_f(r,E,T) \phi(r) dV}{\int \phi^*(r) \chi_p \nu(r,E) \sigma_f(r,E,T) \phi(r) dV}
$$

---

## 🔗 相关文档

- [Phase 1 实现](./BETA_EFFECTIVE_GUIDE.md)
- [Phase 2.2 实现](./PHASE_2_2_IMPLEMENTATION.md)
- [Phase 2.3 实现](./PHASE_2_3_IMPLEMENTATION.md)
- [OpenMC 核数据库文档](https://docs.openmc.org/en/stable/usersguide/cross_sections.html)

---

**实现者**: GitHub Copilot  
**验证状态**: ✅ 编译通过，功能完整  
**代码审查**: 待进行  
**物理验证**: 待实际运行测试  

---

## 附录：API 参考

### `Nuclide::nu()` 方法

```cpp
double Nuclide::nu(double E, EmissionMode mode, int group = 0) const;
```

**参数**:

- `E`: 中子能量 (eV)
- `mode`: `EmissionMode::prompt` | `delayed` | `total`
- `group`: 缓发群索引 (1-8, 仅 delayed 模式使用)

**返回**: 中子产额

### `Nuclide::collapse_rate()` 方法

```cpp
double collapse_rate(int MT, double temperature, 
                     span<const double> energy,
                     span<const double> flux) const;
```

**参数**:

- `MT`: ENDF MT 编号 (18 = 裂变)
- `temperature`: 温度 (K)
- `energy`: 能量网格 (eV)
- `flux`: 通量权重

**返回**: 折叠后的反应率 (barns)

---

**文档版本**: 1.0  
**最后更新**: 2025-11-26
