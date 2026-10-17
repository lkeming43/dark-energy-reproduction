## Algorithms and Methods Used in the Paper

### 1. Data Preprocessing｜数据预处理

- Select DESI $$D_H/r_d$$ data for the DESI-only reconstruction.
  选择DESI的 $$D_H/r_d$$ 数据进行DESI-only重建。
- Select DESI $$D_M/r_d$$ data when combining BAO with supernovae.
  DESI与超新星结合时使用 $$D_M/r_d$$ 数据。
- Remove Pantheon+ data at $$z<0.01$$.
  删除 $$z<0.01$$ 的Pantheon+数据。
- Remove DESY5 supernovae with $$\sigma_m>1$$.
  删除星等误差满足 $$\sigma_m>1$$ 的DESY5超新星。

### 2. Supernova Distance Conversion｜超新星距离转换

Convert apparent magnitude into luminosity distance:

将视星等转换为光度距离：

$$D_L=10^{(m-M-25)/5}\ \mathrm{Mpc}$$

Then convert it into transverse comoving distance:

再转换为横向共动距离：

$$D_M=\frac{D_L}{1+z}$$

Therefore:

因此：

$$\frac{D_M}{r_d}=\frac{10^{(m-M-25)/5}}{(1+z)r_d}$$

### 3. Symbolic Regression｜符号回归

Symbolic regression searches directly for an analytical function:

符号回归直接寻找解析函数：

$$y=f(z)$$

It determines both the mathematical structure and numerical coefficients.

它同时寻找公式结构和数值系数。

### 4. PySR｜高性能符号回归

PySR combines two optimization methods:

PySR结合了两种优化方法：

- Evolutionary algorithms search for the equation structure.
  进化算法搜索公式结构。
- Gradient-based optimization determines numerical constants.
  梯度优化确定公式中的数值常数。

Example:

例如：

$$az\exp(bz)\longrightarrow33.3z\exp(-0.274z)$$

### 5. Operators｜运算符

Unary operators:

一元运算符：

$$\log(x),\quad\sqrt{x},\quad\exp(x),\quad x^2,\quad x^3$$

Binary operators:

二元运算符：

$$+,\quad-,\quad\times,\quad/,\quad x^y$$

### 6. Multi-objective Optimization｜多目标优化

PySR minimizes both prediction error and equation complexity.

PySR同时最小化预测误差和公式复杂度。

The goal is to find an equation that is accurate but still simple and interpretable.

目标是寻找既准确、又简单且可解释的公式。

### 7. Covariance-weighted Loss Function｜协方差加权损失函数

The custom loss function is:

自定义损失函数为：

$$\mathcal{L}=\frac{1}{n}\mathbf{r}^{T}C^{-1}\mathbf{r}$$

where:

其中：

- $$\mathbf{r}$$: residual vector｜残差向量
- $$C$$: covariance matrix｜协方差矩阵
- $$C^{-1}$$: inverse covariance matrix｜逆协方差矩阵
- $$n$$: number of data points｜数据点数量

The residual vector is:

残差向量为：

$$\mathbf{r}=\mathbf{x}*{\mathrm{pred}}-\mathbf{y}*{\mathrm{obs}}$$

This loss accounts for both measurement uncertainties and correlations between observations.

该损失函数同时考虑测量误差和不同观测之间的相关性。

### 8. Model Selection｜模型选择

The paper uses PySR’s `best` selection criterion.

论文采用PySR的 `best` 模型选择标准。

First, equations satisfying

首先保留满足以下条件的公式：

$$\mathcal{L}<1.5\mathcal{L}_{\min}$$

are retained.

Then the equation with the best balance between accuracy and complexity is selected.

然后选择准确度与复杂度平衡最好的公式。

### 9. Monte Carlo Uncertainty Estimation｜蒙特卡洛不确定度估计

Monte Carlo sampling is used to estimate:

使用蒙特卡洛抽样估计：

- $$1\sigma$$ confidence interval｜置信区间
- $$2\sigma$$ confidence interval｜置信区间
- $$3\sigma$$ confidence interval｜置信区间

The method assumes Gaussian reconstruction uncertainties and repeatedly generates possible reconstructed curves.

该方法假设重建误差服从高斯分布，并重复生成可能的重建曲线。

### 10. Differentiation｜求导

For DESI-only:

对于DESI-only：

$$\frac{D_H}{r_d}\longrightarrow\frac{d}{dz}\left(\frac{D_H}{r_d}\right)\longrightarrow w(z)$$

Only the first derivative is required.

只需要一阶导数。

For DESI combined with supernovae:

对于DESI与超新星组合：

$$\frac{D_M}{r_d}\longrightarrow\frac{d}{dz}\left(\frac{D_M}{r_d}\right),\frac{d^2}{dz^2}\left(\frac{D_M}{r_d}\right)\longrightarrow w(z)$$

Both first and second derivatives are required.

需要一阶和二阶导数。

### 11. Reconstruction of $$w(z)$$｜重建 $$w(z)$$

The reconstructed distance function and its derivatives are inserted into the inverted Friedmann equations.

将重建的距离函数及其导数代入反演后的Friedmann方程。

The complete workflow is:

完整流程为：

$$\mathrm{Data}\longrightarrow\mathrm{PySR}\longrightarrow\mathrm{distance\ function}\longrightarrow\mathrm{derivatives}\longrightarrow w(z)$$

### 12. Parameter Sensitivity Test｜参数敏感性检验

The paper tests:

论文检验：

$$\Omega_m=0.25,\ 0.30,\ 0.35$$

and

以及：

$$r_d=90.9,\ 101,\ 111.1\ h^{-1}\mathrm{Mpc}$$

This checks whether the reconstructed $$w(z)$$ depends strongly on the assumed cosmological parameters.

该检验用于判断重建的 $$w(z)$$ 是否强烈依赖假设的宇宙学参数。

### 13. Reference Models｜参考模型

The reconstruction is compared with the $$\Lambda\mathrm{CDM}$$ prediction:

重建结果与 $$\Lambda\mathrm{CDM}$$ 预测比较：

$$w(z)=-1$$

It is also compared with the CPL parameterization:

同时与CPL参数化比较：

$$w(z)=w_0+w_a\frac{z}{1+z}$$

CPL is a reference model, not an algorithm used for reconstruction.

CPL是参考模型，不是论文进行重建的算法。

### Complete Pipeline｜完整流程

$$\mathrm{Data\ cleaning}\longrightarrow\mathrm{distance\ conversion}\longrightarrow\mathrm{PySR}\longrightarrow\mathrm{model\ selection}\longrightarrow\mathrm{Monte\ Carlo}\longrightarrow\mathrm{differentiation}\longrightarrow w(z)$$

$$\mathrm{数据清洗}\longrightarrow\mathrm{距离转换}\longrightarrow\mathrm{PySR}\longrightarrow\mathrm{模型选择}\longrightarrow\mathrm{蒙特卡洛}\longrightarrow\mathrm{求导}\longrightarrow w(z)$$

### Methods Mentioned but Not Used｜提到但未实际使用的方法

- Gaussian Processes (GP)｜高斯过程
- Genetic Algorithms (GA)｜遗传算法
- Eureqa
- gplearn
- AI Feynman
- Barboza-Alcaniz parameterization｜BA参数化

These methods are discussed as background or alternatives but are not used in the main reconstruction.

这些方法只作为背景或替代方案被提及，没有用于论文的主要重建。