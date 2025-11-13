# Cox-Ingersoll-Ross (CIR) Model Estimation and Simulation

### Euler-Maruyama Numerical Method with Monte Carlo Analysis

**Authors:** Group 4  
**Institution:** Kwame Nkrumah University of Science and Technology  
**Course:** Mathematical Finance / Stochastic Processes  
**Year:** 2nd Year Group Project

---

## Abstract

This project implements the **Cox-Ingersoll-Ross (CIR) model** for modeling interest rates, volatility, or other mean-reverting financial processes. We estimate model parameters (κ, θ, σ) using **quasi-maximum likelihood estimation** based on the Euler-Maruyama discretization scheme. The model is then validated through single-path simulation and **Monte Carlo analysis** with 40 independent paths. Our implementation includes numerical stability enhancements (non-negativity enforcement, small epsilon for variance) and comprehensive diagnostics (Feller condition check, convergence analysis, RMSE/MAE metrics).

**Key Result:** Successfully estimated CIR parameters from daily financial data (1998-2024) and demonstrated stable simulation with mean-reverting behavior characteristic of the CIR process.

---

## Repository Contents

```
cir-euler-maruyama-project/
├── README.md                          # This file
├── BESTARIMA_CIR_EulerRevised.R      # Complete R implementation
├── Data.xlsx                            # Input data (daily "High" series)
└── output/
    ├── figures/
    │   ├── 00_series_plot.png         # Original time series
    │   ├── 01_actual_vs_simulated.png # Single-path comparison
    │   └── 02_mc_paths_first10.png    # Monte Carlo simulation (10 paths)
    ├── tables/
    │   ├── data_summary.csv           # Descriptive statistics
    │   ├── estimated_parameters.csv   # κ, θ, σ estimates
    │   ├── single_path_metrics.csv    # RMSE, MAE
    │   ├── param_SEs.csv              # Standard errors (if Hessian available)
    │   └── mc_envelope.csv            # Mean ± 1.96·SD bounds
    └── results/
        └── result_summary.rds         # R object with all results
```

---

## The CIR Model

### Stochastic Differential Equation

The CIR process follows:

**dr(t) = κ(θ − r(t))dt + σ√r(t)dW(t)**

Where:
- **κ (kappa):** Speed of mean reversion (how fast process returns to long-term mean)
- **θ (theta):** Long-term mean level
- **σ (sigma):** Volatility parameter (scales with √r(t))
- **W(t):** Standard Brownian motion

### Key Properties

1. **Mean Reversion:** Process gravitates toward θ over time
2. **Non-negativity:** Volatility √r(t) ensures r(t) ≥ 0 (appropriate for interest rates)
3. **Feller Condition:** If **2κθ ≥ σ²**, process never reaches zero
4. **State-Dependent Volatility:** Higher rates → higher volatility

---

## Implementation Details

### 1. Parameter Estimation (Quasi-MLE)

We maximize the **Euler-based likelihood**:

```r
# For each time step i:
μ_i = κ(θ - r_{i-1})Δt
σ²_i = σ² · r_{i-1} · Δt

# Transition density (approximately Gaussian):
r_i | r_{i-1} ~ N(r_{i-1} + μ_i, σ²_i)

# Negative log-likelihood:
NLL = Σ [0.5·log(2πσ²_i) + (r_i - r_{i-1} - μ_i)² / (2σ²_i)]
```

**Optimization:** L-BFGS-B with positivity constraints on all parameters.

### 2. Euler-Maruyama Discretization

```r
# Numerical scheme for simulation:
r_{i+1} = r_i + κ(θ - r_i)Δt + σ√max(r_i, 0)·√Δt·Z_i

where Z_i ~ N(0,1)
```

**Stability Enhancements:**
- `max(r_i, 0)` prevents negative rates from propagating
- Small epsilon (1e-8) in variance calculations to avoid division by zero
- Explicit non-negativity enforcement after each step

### 3. Monte Carlo Validation

Generate 40 independent paths to assess:
- Variability around mean trajectory
- Long-term behavior (convergence to θ)
- Empirical distribution of terminal values

---

## Quick Start

### Requirements

**R Version:** ≥ 4.0.0

**Required Packages:**
```r
install.packages(c("readxl", "zoo", "ggplot2", "reshape2"))
```

### Running the Analysis

```r
# 1. Set working directory to project folder
setwd("path/to/cir-project")

# 2. Ensure data file exists
# File: g4.xlsx with columns: Date, High

# 3. Run complete analysis
source("BESTARIMA_CIR_EulerRevised.R")

# Expected runtime: ~30 seconds
# Outputs: 3 figures, 5 tables, 1 RData file
```

### Quick Example: Simulating New Paths

```r
# Load estimated parameters
params <- read.csv("output/tables/estimated_parameters.csv")
kappa <- params$kappa
theta <- params$theta
sigma <- params$sigma

# Simulate one new path (1 year, daily steps)
set.seed(456)
n_steps <- 252
dt <- 1/252
r0 <- 20  # starting value

path <- numeric(n_steps + 1)
path[1] <- r0

for(i in 2:(n_steps + 1)) {
  r_prev <- max(path[i-1], 0)
  Z <- rnorm(1)
  path[i] <- max(0, r_prev + kappa*(theta - r_prev)*dt + 
                  sigma*sqrt(r_prev)*sqrt(dt)*Z)
}

# Plot
plot(path, type='l', main="New Simulated CIR Path", 
     ylab="Value", xlab="Time Step")
abline(h=theta, col='red', lty=2)  # long-term mean
```

---

## Results

### Estimated Parameters

Based on daily data from 1998-2024:

| Parameter | Symbol | Estimate | Interpretation |
|-----------|--------|----------|----------------|
| Speed of mean reversion | κ | [your value] | Rate of return to long-term mean |
| Long-term mean | θ | [your value] | Equilibrium level |
| Volatility | σ | [your value] | Proportional to √r(t) |

### Model Fit Quality

| Metric | Value | Interpretation |
|--------|-------|----------------|
| RMSE | [your value] | Average deviation from actual data |
| MAE | [your value] | Mean absolute prediction error |
| Feller Condition | 2κθ - σ² | If > 0: process stays positive |
| Optimization Convergence | 0 | Success (0 = converged) |

### Visual Results

#### Figure 1: Original Time Series (1998-2024)
<img src="output/figures/00_series_plot.png" width="700" alt="Original Series">

**Observation:** Data shows mean-reverting behavior with periods of high/low volatility, making CIR model appropriate.

#### Figure 2: Actual vs. Simulated (Single Path)
<img src="output/figures/01_actual_vs_simulated.png" width="700" alt="Actual vs Simulated">

**Interpretation:**
- **Red:** Actual data with real market shocks and regime changes
- **Cyan:** Simulated path captures general volatility structure and mean reversion
- **Divergence:** Expected - single path won't match exact history, but statistical properties should align

#### Figure 3: Monte Carlo Simulation (10 Sample Paths)
<img src="output/figures/02_mc_paths_first10.png" width="700" alt="Monte Carlo Paths">

**Key Insights:**
- All paths exhibit mean reversion toward θ (not visible here but confirmed in envelope)
- Variability increases when r(t) is higher (state-dependent volatility)
- No negative values (non-negativity enforced)
- Paths diverge but maintain similar statistical properties

---

## Methodology

### Step 1: Data Preparation
```r
# Load daily "High" values from Excel
# Convert to time series with proper date handling
# Calculate summary statistics (mean, variance, range)
```

### Step 2: Parameter Estimation
```r
# Set up negative log-likelihood function
# Use Euler approximation for transition densities
# Optimize with L-BFGS-B (constrained optimization)
# Extract κ, θ, σ and check convergence
```

### Step 3: Model Validation
```r
# 1. Single-path simulation with estimated parameters
# 2. Compare to actual data (RMSE, MAE)
# 3. Check Feller condition: 2κθ ≥ σ²
# 4. Verify optimization convergence (code = 0)
```

### Step 4: Monte Carlo Analysis
```r
# Generate 40 independent paths
# Calculate ensemble mean ± 1.96·SD
# Visualize first 10 paths
# Save confidence envelope for future use
```

---

## Theoretical Background

### Why CIR Model?

**Advantages over Geometric Brownian Motion:**
1. **Mean reversion:** Realistic for interest rates, volatility indices
2. **Non-negativity:** Automatically enforced through √r(t) term
3. **Closed-form solutions:** Bond prices, option formulas available
4. **Empirical support:** Fits short-rate data better than Vasicek model

**Limitations:**
- Assumes constant κ, θ, σ (no time-variation)
- May not capture jumps or regime switches
- Euler discretization introduces bias for large Δt

### Feller Condition Explained

**Condition:** 2κθ ≥ σ²

**Meaning:**
- If satisfied: Process **never** reaches r(t) = 0
- Drift toward θ dominates volatility near zero
- Mathematically: "0 is inaccessible"

**If violated:**
- Process can hit zero with positive probability
- Still non-negative due to √r(t) term
- More frequent boundary encounters

---

## Diagnostics & Validation

### 1. Convergence Check
```r
if(opt$convergence == 0) {
  cat("✓ Optimization converged successfully\n")
} else {
  cat("✗ Warning: convergence code =", opt$convergence, "\n")
}
```

### 2. Parameter Reasonableness
- **κ > 0:** Mean reversion exists (if ≤ 0, no reversion)
- **θ > 0:** Long-term mean is positive
- **σ > 0:** Volatility is present
- **κ not too large:** Avoid over-rapid reversion (numerical instability)

### 3. Simulation Quality
- **Visual inspection:** Paths should "look" mean-reverting
- **No boundary issues:** r(t) stays well above zero most of the time
- **Variance matching:** Simulated variance ~ empirical variance

---

## Applications

### 1. Interest Rate Modeling
- **Short rates:** Fed funds rate, LIBOR
- **Term structure:** CIR++ for multi-factor models
- **Bond pricing:** Closed-form solutions available

### 2. Volatility Modeling
- **Stochastic volatility:** Heston model uses CIR for variance process
- **VIX dynamics:** Mean-reverting volatility index
- **Risk management:** VaR, CVaR calculations

### 3. Commodity Prices
- **Mean reversion:** Oil, gas, agricultural products revert to production costs
- **Storage theory:** Convenience yield follows CIR-like process

---

## Extensions & Future Work

### Immediate Improvements
1. **Exact MLE:** Use true CIR transition density (non-central chi-squared)
2. **Kalman Filter:** For parameter tracking over time
3. **Diagnostic Tests:** Ljung-Box (residual autocorrelation), QQ-plots (normality)

### Advanced Extensions
1. **Jump-Diffusion:** Add Poisson jumps for sudden rate changes
2. **Time-Varying Parameters:** Allow κ(t), θ(t), σ(t) to vary
3. **Multi-Factor Models:** CIR++ with multiple correlated processes
4. **Option Pricing:** Implement bond options, caps, floors

### Research Directions
- Compare Euler vs. Milstein vs. Exact simulation schemes
- Bayesian parameter estimation with MCMC
- Model selection (CIR vs. Vasicek vs. CKLS)

---

## Troubleshooting

### Common Issues

**1. "Column 'High' not found"**
- Ensure Excel file has column named exactly "High"
- Check for extra spaces or different capitalization

**2. "Optimization did not converge"**
- Try different initial values: `init_pars <- c(0.1, median(r), sd(r)*0.5)`
- Increase iterations: `control = list(maxit = 5000)`
- Check for data issues: missing values, extreme outliers

**3. "Simulated paths go negative"**
- Verify `pmax(X_prev, 0)` is applied in simulation loop
- Check if Feller condition violated (expect boundary encounters)
- Consider smaller Δt for better discretization

**4. "Poor fit (high RMSE)"**
- CIR assumes constant parameters - data may have regime changes
- Try estimating on subperiods separately
- Consider richer models (jump-diffusion, time-varying parameters)

---

## File Descriptions

| File | Description | Size |
|------|-------------|------|
| `BESTARIMA_CIR_EulerRevised.R` | Complete implementation with comments | ~200 lines |
| `g4.xlsx` | Daily price data (Date, High columns) | Variable |
| `00_series_plot.png` | Original time series visualization | ~50 KB |
| `01_actual_vs_simulated.png` | Model fit comparison | ~100 KB |
| `02_mc_paths_first10.png` | Monte Carlo sample paths | ~80 KB |
| `estimated_parameters.csv` | κ, θ, σ values | <1 KB |
| `single_path_metrics.csv` | RMSE, MAE | <1 KB |
| `mc_envelope.csv` | Mean ± confidence bands | ~10 KB |

---

## Citation

```
Group 4 (2025). Cox-Ingersoll-Ross Model Estimation and Simulation 
using Euler-Maruyama Method. Department of Statistics, 
Kwame Nkrumah University of Science and Technology.
```

**BibTeX:**
```bibtex
@misc{group4_2025_cir,
  title={Cox-Ingersoll-Ross Model Estimation and Simulation using Euler-Maruyama Method},
  author={Group 4},
  year={2025},
  institution={Kwame Nkrumah University of Science and Technology},
  howpublished={GitHub Repository}
}
```

---

## References

### Core Papers
1. **Cox, J. C., Ingersoll, J. E., & Ross, S. A. (1985).** "A Theory of the Term Structure of Interest Rates." *Econometrica*, 53(2), 385-407.
2. **Kloeden, P. E., & Platen, E. (1992).** *Numerical Solution of Stochastic Differential Equations.* Springer.

### Implementation Guides
3. **Aït-Sahalia, Y. (2002).** "Maximum Likelihood Estimation of Discretely Sampled Diffusions: A Closed-Form Approximation Approach." *Econometrica*, 70(1), 223-262.
4. **Glasserman, P. (2004).** *Monte Carlo Methods in Financial Engineering.* Springer.

---

## License

This project is available for academic and educational use. Please cite appropriately if used in research.

---

## Contact

**Bentum Welson (Group leader)**  
Department of Statistics  
Kwame Nkrumah University of Science and Technology, Ghana

For questions, please open an issue in this repository or contact group members.

---

**Keywords:** CIR Model · Euler-Maruyama · Stochastic Differential Equations · Interest Rate Modeling · Monte Carlo Simulation · Maximum Likelihood Estimation · Mean Reversion · Mathematical Finance · R Programming
