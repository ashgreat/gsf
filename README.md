# GSF: Generalized Stochastic Frontier Analysis

This R package implements the Generalized Stochastic Frontier (GSF) model with input-specific inefficiencies as described in Kumbhakar and Tsionas (2021). It uses `rstan` for Bayesian estimation.

## Installation

You can install the package from GitHub using the `remotes` package:

```R
# install.packages("remotes")
remotes::install_github("yourusername/GSF")
```

## Usage

Here is a simple example using simulated data:

```R
library(GSF)

# 1. Simulate data
set.seed(123)
sim_data <- simulate_gsf_data(N = 100, J = 2, M = 2)

# 2. Fit the model
fit <- gsf(
  y = sim_data$y,
  X = sim_data$X,
  Z = sim_data$Z,
  iter = 2000,
  chains = 4
)

# 3. Print summary
print(fit, pars = c("alpha", "sigma_sq", "sigma_u_sq"))

# 4. Extract inefficiencies
ineff <- extract_inefficiencies(fit)

# Input slacks (negative values, closer to 0 is more efficient)
head(ineff$input_slacks)

# Output inefficiency (negative values)
head(ineff$output_inefficiency)
```

## Model Details

The model estimates:
- **Input Slacks ($\theta$)**: The inefficiency associated with each input.
- **Output Inefficiency ($u_0$)**: The overall technical inefficiency.

The production function is specified as a translog function of effective inputs ($x + \theta$).

## References

Kumbhakar, S. C., & Tsionas, M. G. (2021). Dissections of input and output efficiency: A generalized stochastic frontier model. *International Journal of Production Economics*, 232, 107940.
