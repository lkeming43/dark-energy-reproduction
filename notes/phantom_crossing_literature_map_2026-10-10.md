# “不同来源的 \(w=-1\) 穿越 + 假穿越校准”文献地图

检索日期：2026-10-10（Asia/Singapore）

## 1. 结论先行

截至本次检索，没有发现一篇论文同时完成以下四件事：

1. 用 DESI DR2 BAO 与 Pantheon+/DESY5 的实际红移位置和协方差生成模拟数据；
2. 用 PySR 或同一套符号回归流程重建距离、\(H(z)\) 和 \(w(z)\)；
3. 定量测量 \(\Lambda\)CDM 下的“假穿越率”和“假振荡率”，以及注入真实穿越后的检出率和区间覆盖率；
4. 在统一数据、噪声和分析流程下，比较 CPL、相互作用暗部门、质量变化中微子和修正引力产生的相似有效 \(w(z)\)，研究仅靠背景距离能否识别物理来源。

但是，这四部分分别都有很接近的前人工作。因此合适的定位是：**把已有的符号回归、模拟校准和暗简并研究接到同一条 DESI DR2 + SN 流水线上**，而不是声称首次研究 phantom crossing。

## 2. 与拟议项目最接近的论文

| 优先级 | 论文 | 已经做了什么 | 与拟议项目之间的空缺 |
|---|---|---|---|
| 必读 | [Sousa-Neto et al., 2025, arXiv:2502.10506](https://arxiv.org/abs/2502.10506) | 用 PySR 从 DESI DR2、Pantheon+、DESY5 重建 \(w(z)\)；在所考察范围内没有超过 \(3\sigma\) 的偏离 | 没有用大量已知真值的 mock 计算假穿越率、检出率和覆盖率，也不识别物理来源 |
| 必读 | [Borghetto et al., 2026, arXiv:2606.17951](https://arxiv.org/abs/2606.17951) | 先做贝叶斯样条重建，再用 exhaustive symbolic regression 找到简单穿越形式 \(w(a)=w_0/\sqrt a\) | 重点是发现低复杂度参数化；没有针对 PySR 的 \(\Lambda\)CDM 假穿越校准，也没有多种物理来源的混淆实验 |
| 必读 | [Arjona, Melchiorri & Nesseris, 2021, arXiv:2107.04343](https://arxiv.org/abs/2107.04343) | 为 \(\Lambda\)CDM、\(w\)CDM 和两种 \(f(R)\) 模型生成 LSST 风格增长率 mock，并用遗传算法重建零检验 | 方法论最接近“已知真值 → 加噪声 → ML 重建 → 判断能否区分”；使用增长数据和 GA，不是 DESI+SN、PySR 或 crossing 统计 |
| 必读 | [Jiang, Hong & Zhang, 2026, arXiv:2605.13427](https://arxiv.org/abs/2605.13427) | 用 \(\Lambda\)CDM 与轻微 \(w_0w_a\) 演化的 mock 检查 GP 重建，并比较把先验放在 \(H\)、\(f\)、\(\log f\) 上的影响 | 已做注入恢复，但没有符号回归，也没有直接报告 crossing 假阳性率 |
| 必读 | [Nandi, Ghoshal & Mota, 2026, arXiv:2605.27357](https://arxiv.org/abs/2605.27357) | 为六种背景模型生成未来标准汽笛 mock，用 GP 重建并定量比较模型可分辨性；结论是背景距离通常不能决定性地区分模型 | 最接近“不同物理来源是否可识别”；数据是未来引力波标准汽笛，方法是 GP，没有 DESI+SN/PySR 的假穿越校准 |
| 必读 | [von Marttens et al., 2019, arXiv:1911.02618](https://arxiv.org/abs/1911.02618) | 明确推导动力学暗能量与相互作用暗能量在背景和扰动层面的暗简并条件 | 提供理论基线；没有符号回归或 mock 混淆矩阵 |
| 高 | [Li & Wang, 2025, arXiv:2506.22953](https://arxiv.org/abs/2506.22953) | DESI DR2+CMB+三套 SN，用分箱和多项式插值重建 \(w(z)\) 与暗能量密度 | 可作为 PySR 的公平方法对照；未做重复 mock 的假穿越校准 |
| 高 | [Mukherjee & Sen, 2024, arXiv:2405.19178](https://arxiv.org/abs/2405.19178) | GP 重建 DESI DR1+不同 SN；不同 SN 组合给出明显不同结论 | 直接说明数据组合可制造或消除穿越信号；没有统一 mock 真值测试 |
| 高 | [Calderon et al., 2024, arXiv:2405.04216](https://arxiv.org/abs/2405.04216) | 用 Crossing Statistics 对 DESI DR1+CMB+三套 SN 做模型无关重建 | 是另一种 reconstruction baseline；没有 PySR 算法假阳性校准 |
| 高 | [Koksbang & Heinesen, 2026, arXiv:2604.05822](https://arxiv.org/abs/2604.05822) | 用 bootstrap symbolic regression 重建距离、\(H\) 及其导数，检验广义 FLRW 一致性 | 证明符号回归与 bootstrap 可以结合；目标不是 \(w=-1\) 穿越和物理来源识别 |

## 3. 符号回归、遗传算法与机器学习重建

1. [Nesseris & Garcia-Bellido, “A new perspective on Dark Energy modeling via Genetic Algorithms,” arXiv:1205.0364](https://arxiv.org/abs/1205.0364)  
   遗传算法重建距离、增长率和 \(w(z)\)，并用路径积分、Fisher 和 bootstrap 比较误差估计。它是当前 PySR 工作的重要方法学祖先。

2. [Gangopadhyay, Sami & Sharma, “Phantom dark energy as a natural selection of evolutionary processes à la genetic algorithm…,” arXiv:2303.07301](https://arxiv.org/abs/2303.07301)  
   用遗传算法结合低红移背景和 RSD 数据，得到偏 phantom 的重建。适合作为“算法可以选出 phantom 形状”的先例。

3. [Sousa-Neto et al., “Evidence for dynamical dark energy from DESI-DR2 and SN data? A symbolic regression analysis,” arXiv:2502.10506](https://arxiv.org/abs/2502.10506)  
   当前复现项目的直接对象。

4. [Borghetto et al., “Data-Driven Discovery of a Simple Phantom-Crossing Dark Energy Parametrization,” arXiv:2606.17951](https://arxiv.org/abs/2606.17951)  
   目前与“symbolic regression + phantom crossing”最直接的另一篇论文。

5. [Arjona, Melchiorri & Nesseris, “Testing the \(\Lambda\)CDM paradigm with growth rate data and machine learning,” arXiv:2107.04343](https://arxiv.org/abs/2107.04343)  
   用多种宇宙模型的 mock 验证遗传算法零检验。

6. [Koksbang & Heinesen, “Model-independent constraints on generalized FLRW consistency relations with bootstrap-based symbolic regression,” arXiv:2604.05822](https://arxiv.org/abs/2604.05822)  
   对符号回归的不确定度和导数稳定性特别有参考价值。

## 4. 非参数重建、模拟注入与“如何避免把噪声当结构”

1. [Clarkson & Zunckel, “Direct reconstruction of dark energy,” arXiv:1002.5004](https://arxiv.org/abs/1002.5004)  
   用 PCA 和信息准则识别真实结构并平滑噪声，属于“假振荡”问题的早期直接方法论文。

2. [Holsclaw et al., “Nonparametric Reconstruction of the Dark Energy Equation of State,” arXiv:1009.5443](https://arxiv.org/abs/1009.5443)  
   在模拟超新星数据上验证 GP 重建非平凡 \(w(z)\) 和误差区间。

3. [Holsclaw et al., “Nonparametric Dark Energy Reconstruction from Supernova Data,” arXiv:1011.3079](https://arxiv.org/abs/1011.3079)  
   将 GP+MCMC 用于真实超新星数据。

4. [Holsclaw et al., “Nonparametric Reconstruction … from Diverse Data Sets,” arXiv:1104.2041](https://arxiv.org/abs/1104.2041)  
   扩展到多种宇宙学探针。

5. [Seikel, Clarkson & Smith, “Reconstruction of dark energy and expansion dynamics using Gaussian processes,” arXiv:1204.2832](https://arxiv.org/abs/1204.2832)  
   GP 重建及导数误差传播的经典论文，并在预报数据上测试可恢复精度。

6. [Gerardi, Martinelli & Silvestri, “Reconstruction … impact of theoretical priors,” arXiv:1902.09423](https://arxiv.org/abs/1902.09423)  
   展示 quintessence/Horndeski 等理论先验如何改变 \(w(z)\) 重建。

7. [Lin, Li & Tang, “Non-parametric reconstruction … from Pantheon,” arXiv:1905.11593](https://arxiv.org/abs/1905.11593)  
   强调 \(w(z)\) 依赖距离的二阶导数，因此误差会被放大，并对 \(H_0\)、\(\Omega_m\) 先验敏感。

8. [Zhang et al., “Kernel Selection for Gaussian Process in Cosmology,” arXiv:2304.03911](https://arxiv.org/abs/2304.03911)  
   说明 GP 核函数本身也是分析选择，适合放进方法敏感性章节。

9. [Escamilla et al., “Model-independent reconstruction of the Interacting Dark Energy Kernel,” arXiv:2305.16290](https://arxiv.org/abs/2305.16290)  
   分箱和 GP 重建得到约 \(1\sigma\) 的振荡迹象，但贝叶斯证据惩罚额外自由度；这是“曲线有波纹不等于新物理”的好例子。

10. [You, Wang & Yang, “Dynamical Dark Energy Implies a Coupled Dark Sector…,” arXiv:2504.00985](https://arxiv.org/abs/2504.00985)  
    GP 重建 DESI DR2 背景并进一步推断暗部门耦合。

11. [Ye et al., “Late-time reconstruction of non-minimally coupled gravity with a smoothness prior,” arXiv:2605.12415](https://arxiv.org/abs/2605.12415)  
    用相关平滑先验约束高自由度函数，是设计“复杂度/平滑度正则化”对照的关键参考。

12. [Jiang, Hong & Zhang, arXiv:2605.13427](https://arxiv.org/abs/2605.13427)  
    直接包含 \(\Lambda\)CDM 和轻微演化的 mock 注入恢复测试。

13. [Nandi, Ghoshal & Mota, arXiv:2605.27357](https://arxiv.org/abs/2605.27357)  
    直接包含六种物理背景的 mock 和模型可分辨性度量。

## 5. 为什么同一个“穿越”可能来自不同物理

### 5.1 暗简并与解释框架

1. [Kunz, “The dark degeneracy,” arXiv:astro-ph/0702615](https://arxiv.org/abs/astro-ph/0702615)  
   核心结论：引力首先测量总能动张量；仅靠几何膨胀史不能唯一拆分暗物质和暗能量。

2. [von Marttens et al., “Dark degeneracy I,” arXiv:1911.02618](https://arxiv.org/abs/1911.02618)  
   系统推导动力学暗能量和相互作用暗能量的等价条件，并讨论扰动信息如何部分破除简并。

3. [Mishra, “Effective Phantom Dark Energy: What Cosmological Reconstruction Does and Does Not Imply,” arXiv:2605.27301](https://arxiv.org/abs/2605.27301)  
   专门澄清：重建出的有效 \(w<-1\) 并不自动等于基本 phantom 场、鬼不稳定或未来 Big Rip。

4. [Kunz & Sapone, “Dark Energy versus Modified Gravity,” arXiv:astro-ph/0612452](https://arxiv.org/abs/astro-ph/0612452)  
   展示广义暗能量流体可以模拟修正引力在背景和扰动上的行为，说明识别问题本身很难。

### 5.2 质量变化中微子

1. [Fardon, Nelson & Weiner, “Dark Energy from Mass Varying Neutrinos,” arXiv:astro-ph/0309800](https://arxiv.org/abs/astro-ph/0309800)  
   MaVaN 思路的奠基论文。

2. [Brookfield et al., “Cosmology of Mass-Varying Neutrinos Driven by Quintessence,” arXiv:astro-ph/0512367](https://arxiv.org/abs/astro-ph/0512367)  
   已明确指出：中微子与 quintessence 耦合可使超新星推断的表观 \(w<-1\)。

3. [Rodrigues et al., “Effective Neutrino Mass and Apparent Phantom Crossing,” arXiv:2610.10195](https://arxiv.org/abs/2610.10195)  
   直接用 CMB+DESI DR2+SN 研究质量变化中微子造成的表观穿越；是物理来源模拟最相关的近期论文。

### 5.3 相互作用暗部门

1. [Khoury, Lin & Trodden, “Apparent \(w<-1\) and a Lower \(S_8\) from Dark Axion and Dark Baryons Interactions,” arXiv:2503.16415](https://arxiv.org/abs/2503.16415)  
   变化的暗物质质量产生表观 phantom crossing，并影响结构增长。

2. [Guedezounme, Dinda & Maartens, “Phantom crossing or dark interaction?,” arXiv:2507.18274](https://arxiv.org/abs/2507.18274)  
   直接比较“内禀非 phantom 暗能量 + 相互作用”与有效 crossing 的解释。

3. [Antusch, King & Wang, “Coupled Dark Energy and Dark Matter for DESI,” arXiv:2604.08449](https://arxiv.org/abs/2604.08449)  
   区分标量场自身的 \(w_\phi\) 与错误假设暗部门不相互作用时推断出的 \(w_{\rm eff}\)。

4. [Khoury, Lin & Trodden, “Cosmological Evidence for Dark Axion-Dark Baryon Interactions from Apparent Phantom Crossing,” arXiv:2607.16191](https://arxiv.org/abs/2607.16191)  
   将具体相互作用模型写进 Boltzmann 代码并与 CMB+DESI DR2+SN 拟合。

### 5.4 修正引力

1. [Bamba, “Equation of State for Dark Energy in Modified Gravity Theories,” arXiv:1202.4317](https://arxiv.org/abs/1202.4317)  
   综述多类修正引力如何产生有效 phantom divide crossing。

2. [Chudaykin & Kunz, “Modified gravity interpretation of the evolving dark energy in light of DESI data,” arXiv:2407.02558](https://arxiv.org/abs/2407.02558)  
   用 Horndeski 类物理解释 DESI 偏好的演化暗能量，并指出特定 DESI LRG 点的影响。

3. [Efstratiou, Paraskevas & Perivolaropoulos, “Addressing the DESI DR2 Phantom-Crossing Anomaly … with Reconstructed Scalar-Tensor Gravity,” arXiv:2511.04610](https://arxiv.org/abs/2511.04610)  
   用标量张量引力重建同时讨论 crossing 与 \(H_0\)。

4. [Ye et al., arXiv:2605.12415](https://arxiv.org/abs/2605.12415)  
   非最小耦合引力的非参数重建与平滑先验。

5. [Hemmatian et al., “Gravitational-Susceptibility Dip Behind DESI's results…,” arXiv:2609.26517](https://arxiv.org/abs/2609.26517)  
   将类似的背景行为映射到修改后的引力响应。

### 5.5 真正跨越 \(w=-1\) 的场论与扰动处理

1. [Hu, “Crossing the Phantom Divide: Dark Energy Internal Degrees of Freedom,” arXiv:astro-ph/0410680](https://arxiv.org/abs/astro-ph/0410680)  
   说明单一有效流体直接跨越时的扰动困难，以及额外内部自由度如何解决。

2. [Kunz & Sapone, “Crossing the Phantom Divide,” arXiv:astro-ph/0609040](https://arxiv.org/abs/astro-ph/0609040)  
   系统讨论穿越附近压力扰动的有限性以及两场 quintom 示例。

3. [Di Valentino, Mukherjee & Sen, “Dark Energy with Phantom Crossing and the \(H_0\) tension,” arXiv:2005.12587](https://arxiv.org/abs/2005.12587)  
   参数化 crossing 与 \(H_0\) 张力的观测应用。

## 6. 统计解释和稳健性

1. [Zhao et al., “Dynamical dark energy in light of the latest observations,” arXiv:1701.08165](https://arxiv.org/abs/1701.08165)  
   很清楚地区分“最佳拟合改善达到若干 \(\sigma\)”与“贝叶斯证据是否支持更复杂模型”。

2. [Patel, Chakraborty & Amendola, “The prior dependence of the DESI results,” arXiv:2407.06586](https://arxiv.org/abs/2407.06586)  
   改变 \(w_0,w_a\) 先验边界可改变模型比较结论。

3. [Silva & Nunes, “Testing signatures of phantom crossing through full-shape galaxy clustering analysis,” arXiv:2507.13989](https://arxiv.org/abs/2507.13989)  
   加入 BOSS 全形状聚类后报告更强 crossing 信号；可用来说明增长/全形状信息可能帮助破简并，但结论仍依赖参数化。

4. [DESI Collaboration et al., “Constraints on Physics-Focused Aspects of Dark Energy using DESI DR1 BAO Data,” arXiv:2405.13588](https://arxiv.org/abs/2405.13588)  
   对多类有物理动机的暗能量模型做统一数据比较，是物理模型选择的基准。

## 7. 对当前项目最准确的创新表述

建议不要写“首次研究 \(w=-1\) 穿越的不同物理来源”，因为暗简并、相互作用、MaVaN 和修正引力文献已经很多。

可以写成下面这种更可检验的表述：

> 在与 DESI DR2 和超新星观测相匹配的数据协方差下，首次系统校准符号回归重建中 phantom crossing 和振荡特征的假阳性率、检出率及置信区间覆盖率，并量化仅用背景距离时，不同物理机制产生的有效 \(w(z)\) 在同一 PySR 流水线中的可辨识程度。

这句话是否最终能使用，还需要在正式写作前做一次 ADS/INSPIRE 引文追踪和更新检索；当前检索支持它是一个明确的文献空缺，但互联网检索无法证明“绝对没有任何遗漏”。

## 8. 最小可行分析与完整扩展的边界

### 最小可行分析（应先完成）

1. 只生成两类真值：\(\Lambda\)CDM（绝不穿越）和一个会穿越的 CPL。
2. 使用真实 DESI DR2+SN 的红移采样和完整协方差，分别生成至少数百个 mock。
3. 对每个 mock 运行完全相同的 PySR 设置和固定的表达式选择规则。
4. 预先定义 crossing 判据，例如：在规定红移区间内，置信带是否支持从 \(w>-1\) 到 \(w<-1\) 的连续转换。
5. 报告：假阳性率、真阳性率、穿越红移偏差、置信区间覆盖率、对随机种子和复杂度惩罚的敏感性。

这一步回答“PySR 画出的穿越有多可信”，无需先实现复杂微观模型。

### 第二阶段

加入相互作用暗部门或质量变化中微子的一种代表模型，以及一种修正引力模型。先令它们在背景距离上尽量接近同一个 crossing CPL，再比较 PySR 输出。如果背景数据无法区分，这是有意义的负结果；随后加入 RSD/弱透镜等增长信息，检查可辨识度是否提高。

