# gsf: Generalized Stochastic Frontier Model

An R package implementing the Generalized Stochastic Frontier (GSF) model from Kumbhakar and Tsionas (2021) "Dissections of input and output efficiency: A generalized stochastic frontier model." The package now wraps the Bayesian backend with `brms`/Stan to provide familiar posterior-analysis tooling and compatibility with the broader `brms` ecosystem.

## Overview

The `gsf` package jointly estimates:

1. **Output technical inefficiency** (u₀): Percentage reduction in output due to managerial inefficiency
2. **Input-specific inefficiencies/slacks** (θⱼ): Percentage over-use of each input

This provides more detailed efficiency analysis than standard stochastic frontier models by identifying which inputs are being used inefficiently.

## Features

- Translog production function with flexible substitution patterns
- Joint estimation of output and input inefficiencies
- Exogenous determinants of input slacks (z variables)
- Full Bayesian inference via Markov Chain Monte Carlo (Stan) with a lightweight `brmsfit` interface for downstream diagnostics
- Comprehensive diagnostic and visualization tools
- Functions for efficiency extraction and model comparison

## Installation

### Prerequisites

You need to have Stan installed. Install RStan first:

```r
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
```

### Install gsf

From source:

```r
# Install dependencies
install.packages(c("rstan", "brms", "ggplot2", "tidyr", "loo"))

# Or install via remotes from GitHub
remotes::install_github("ashgreat/gsf")
```

## Quick Start

```r
library(gsf)

# Simulate UK manufacturing-like data
uk_data <- simulate_uk_manufacturing(
  n_firms = 100,
  avg_periods = 9,
  seed = 42
)

# Fit the GSF model
fit <- gsf_fit(
  formula = log(Y) ~ log(L) + log(K),
  data = uk_data,
  z_vars = c("CR", "MKSH", "FP", "IMP", "TREND"),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

# View results
summary(fit)

# Efficiency summary (Table 2 from paper)
efficiency_summary(fit)

# Extract inefficiency measures
ti <- extract_technical_inefficiency(fit)
slacks <- extract_input_slacks(fit)

# Visualizations
plot(fit, type = "efficiency")
plot(fit, type = "slacks")
plot_marginal_effects(fit)
```

## Model Specification

The GSF model uses a translog production function:

```
y_i = α₀ + δ₀'z̃_i + α'(x_i + θ_i) + ½(x_i + θ_i)'B(x_i + θ_i) + u₀_i + ε₀_i
```

Where:
- `y_i = ln Y_i` is log output
- `x_i = ln X_i` is the vector of log inputs
- `θ_i ≤ 0` are input-specific slacks (negative = over-use)
- `u₀_i ≤ 0` is output technical inefficiency
- `ε₀_i ~ N(0, σ²)` is random noise

Input slacks follow:
```
θ_i = Δz_i + u_i, u_i ~ N_J(0, Ω), θ_i ≤ 0
```

## Key Functions

| Function | Description |
|----------|-------------|
| `gsf_fit()` | Main model fitting function |
| `extract_technical_inefficiency()` | Extract output technical inefficiency |
| `extract_input_slacks()` | Extract input-specific slacks |
| `extract_elasticities()` | Extract input elasticities |
| `extract_returns_to_scale()` | Extract returns to scale |
| `extract_output_loss()` | Extract output loss estimates |
| `efficiency_summary()` | Summary table of efficiency measures |
| `extract_delta_coefficients()` | Slack determinant coefficients |
| `plot.gsf_fit()` | Various diagnostic plots |
| `plot_marginal_effects()` | Effects of z variables on slacks |
| `waic()` | Model comparison via WAIC |
| `loo_cv()` | Leave-one-out cross-validation |

## Example Output

Based on the UK manufacturing example from the paper:

```
Efficiency Measures (Sample Averages):
                          measure    mean     sd
Technical Inefficiency (u0)    0.024  0.018
Slack (log(L))                 0.023  0.038
Slack (log(K))                 0.107  0.029
Elasticity (log(L))            0.519  0.021
Elasticity (log(K))            0.414  0.042
Returns to Scale               0.933  0.131
```

## Citation

If you use this package, please cite:

```
Kumbhakar, S.C. and Tsionas, M.G. (2021). "Dissections of input and output
efficiency: A generalized stochastic frontier model." International Journal
of Production Economics, 232, 107940.
```

## License

MIT License

## Related Work

This package implements the methodology from:

- Kumbhakar and Tsionas (2021) - Main GSF model
- Kumbhakar (1988) - Input-specific technical inefficiency
- Kopp (1981) - Graphical concept of ISTI
- Koop et al. (2000) - Effective input idea

## Contact

For issues and feature requests, please open an issue on GitHub.
