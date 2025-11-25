library(rstan)
library(MASS)

# Load the package functions (source them for now as package is not installed)
source("R/gsf_model.R")
source("R/simulation.R")
source("R/utils.R")

# Simulate data
set.seed(123)
sim_data <- simulate_gsf_data(N = 50, J = 2, M = 2)

# Fit model
# Use short chains for testing
fit <- gsf(
    y = sim_data$y,
    X = sim_data$X,
    Z = sim_data$Z,
    iter = 1000,
    chains = 2,
    cores = 1
)

# Check results
print(fit, pars = c("alpha0", "alpha", "sigma_sq", "sigma_u_sq"))

# Extract inefficiencies
ineff <- extract_inefficiencies(fit)
head(ineff$input_slacks)
head(ineff$output_inefficiency)

# Compare with true latent values
plot(as.vector(sim_data$latent$theta), as.vector(ineff$input_slacks),
    main = "True vs Estimated Input Slacks", xlab = "True", ylab = "Estimated"
)
abline(0, 1, col = "red")
