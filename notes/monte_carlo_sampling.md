# Monte Carlo Sampling Notes

## 1. Purpose

Monte Carlo sampling is used to propagate the uncertainty of the original observations to the final fitted function or reconstructed result.

Suppose the dataset contains six measurements:

$$\mathbf{z}=(z_1,z_2,\ldots,z_6)$$

$$\mathbf{y}=(y_1,y_2,\ldots,y_6)$$

The redshift positions $\mathbf{z}$ normally remain fixed. Monte Carlo sampling changes the measured values $\mathbf{y}$ according to their uncertainties.

Each Monte Carlo realization still contains only six data points. Repeating the sampling many times does not create new independent observations; it only creates many possible versions of the original six measurements.

## 2. Sampling Methods

### 2.1 Independent-error sampling

If only the standard deviation $\sigma_i$ of each measurement is known, and the six measurements are assumed to be independent, the $i$th simulated value in realization $k$ is

$$y_i^{(k)}=y_i+\sigma_i\epsilon_i^{(k)}$$

where

$$\epsilon_i^{(k)}\sim\mathcal{N}(0,1)$$

Equivalently,

$$y_i^{(k)}\sim\mathcal{N}(y_i,\sigma_i^2)$$

Python example:

```python
import numpy as np

rng = np.random.default_rng(42)

y_mc = rng.normal(
    loc=y,
    scale=sigma,
    size=(n_realizations, 6),
)
```

The resulting array has shape `(n_realizations, 6)`. Each row is one possible realization of the six measurements.

### 2.2 Covariance-aware sampling

If the measurements are correlated and a covariance matrix $C$ is available, the six values must be sampled jointly:

$$\mathbf{y}^{(k)}\sim\mathcal{N}(\mathbf{y},C)$$

Here, $C$ is a $6\times6$ matrix. Its diagonal elements contain the variance of each measurement:

$$C_{ii}=\sigma_i^2$$

Its off-diagonal elements describe how different measurements vary together:

$$C_{ij}=\operatorname{Cov}(y_i,y_j)$$

Python example:

```python
y_mc = rng.multivariate_normal(
    mean=y,
    cov=C,
    size=n_realizations,
)
```

The same sampling can be written using a Cholesky decomposition:

$$C=LL^T$$

$$\mathbf{y}^{(k)}=\mathbf{y}+L\boldsymbol{\epsilon}^{(k)}$$

where

$$\boldsymbol{\epsilon}^{(k)}\sim\mathcal{N}(\mathbf{0},I)$$

Covariance-aware sampling should be used when a reliable covariance matrix is provided. Sampling every point independently would otherwise discard known correlations.

## 3. Overall Workflow

1. Prepare the original measurement locations $\mathbf{z}$ and measured values $\mathbf{y}$.
2. Obtain the measurement uncertainties $\boldsymbol{\sigma}$ or the full covariance matrix $C$.
3. Choose a probability distribution for the observational errors.
4. Choose the number of Monte Carlo realizations $N_{\mathrm{MC}}$.
5. Choose and record a random seed.
6. Generate one simulated six-point dataset $\mathbf{y}^{(k)}$.
7. Fit the model or run PySR using $(\mathbf{z},\mathbf{y}^{(k)})$.
8. Evaluate the fitted function on a common redshift grid.
9. Repeat Steps 6–8 for all $N_{\mathrm{MC}}$ realizations.
10. At each redshift, summarize the distribution of reconstructed values using the median, 68% interval, or 95% interval.
11. Check that the result is stable when $N_{\mathrm{MC}}$ or the random seed is changed.

The complete process is

$$\text{observations}\rightarrow\text{error model}\rightarrow\text{random realizations}\rightarrow\text{repeated fitting}\rightarrow\text{distribution of reconstructed results}$$

## 4. New Assumptions and User-chosen Parameters

Monte Carlo sampling introduces choices that are not additional observational information. They must be stated explicitly because they affect the result.

### 4.1 Probability distribution

The probability distribution specifies how simulated measurements are generated around the observed values.

A common assumption is a Gaussian distribution:

$$\mathbf{y}^{(k)}\sim\mathcal{N}(\boldsymbol{\mu},C)$$

This assumes that the observational errors can be adequately represented by a multivariate normal distribution. If the true errors are asymmetric, bounded, or strongly non-Gaussian, a Gaussian model may be inappropriate.

### 4.2 Distribution center

The center $\boldsymbol{\mu}$ determines the values around which the simulations fluctuate. The usual choice is

$$\boldsymbol{\mu}=\mathbf{y}_{\mathrm{observed}}$$

This treats the published measurements as the centers of their possible repeated observations.

Using a theoretical prediction or a fitted curve as the center would answer a different question and must be documented separately.

### 4.3 Uncertainty scale

For independent sampling, the uncertainty scale is given by the standard deviations:

$$\boldsymbol{\sigma}=(\sigma_1,\sigma_2,\ldots,\sigma_6)$$

For correlated sampling, the uncertainty structure is specified by the covariance matrix $C$.

These values should normally come from the original survey or data release rather than being chosen arbitrarily.

### 4.4 Number of realizations

$N_{\mathrm{MC}}$ is the number of simulated datasets. A larger value produces more stable numerical summaries but does not add new observational information.

The Monte Carlo numerical error generally decreases approximately as

$$\mathrm{error}\propto\frac{1}{\sqrt{N_{\mathrm{MC}}}}$$

Therefore, reducing the Monte Carlo error by a factor of two usually requires about four times as many realizations.

### 4.5 Random seed

The random seed determines the pseudorandom number sequence:

```python
rng = np.random.default_rng(42)
```

It is used for reproducibility. The seed does not make the method more accurate and should not materially change the final result when enough realizations are used.

### 4.6 Treatment of the redshift values

The usual assumption is that the redshifts are fixed:

$$z_i^{(k)}=z_i$$

Only the measured quantities are sampled. If redshift uncertainty is significant, a separate error model for $z_i$ is required.

### 4.7 Model-refitting rule

The analysis must specify what is repeated for each realization:

- **Fixed functional structure:** keep the equation form fixed and refit only its numerical parameters.
- **Full PySR retraining:** allow PySR to search for a new equation structure for every realization.

The first approach mainly propagates measurement and parameter uncertainty. The second also includes variation in the discovered equation structure, but is more computationally expensive and more difficult to summarize.

## 5. Essential Interpretation

Monte Carlo sampling does not turn six measurements into a larger observational dataset. It asks:

> If the same six quantities were measured again under the stated error model, how much could the fitted or reconstructed result change?

The validity of the output depends on the validity of the selected probability distribution, center, uncertainty information, correlations, and refitting procedure.
