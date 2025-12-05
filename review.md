# Analysis of GSF Model Implementation vs. Kumbhakar & Tsionas (2021) Paper

## Summary

I found **several critical issues** with the implementation that need to be corrected to match the paper exactly.

---

## Critical Issues

### 1. CRITICAL: Truncation Constraint on Input Slacks θ is INCORRECTLY Implemented

**Paper specification (Equation 3b, 3c, page 4):**

```
θ_i ≡ ln Θ_i = Δz_i + u_i,  θ_i ≤ 0_{(J×1)}
```

with **Assumption 4b**: `θ_i|z_i ~ N_J(Δz_i, Ω), θ_i ≤ 0_{(J×1)}`

This is a **truncated multivariate normal distribution** where the truncation is enforced at θ ≤ 0.

**Current implementation in `gsf_model.stan` (lines 88-91):**

```stan
for (n in 1:N) {
    vector[J] mu_theta = Delta * Z_slack[n]';
    theta[n]' ~ multi_normal_cholesky(mu_theta, L_Sigma);
}
```

**Problem:** The code uses an **ordinary** multivariate normal distribution, NOT a truncated one. The constraint `matrix<upper=0>[N, J] theta;` only constrains the parameter space but does NOT properly account for the **normalizing constant** of the truncated distribution.

According to the paper (page 5), the pdf of θ_i is:

```
p(θ_i|z_i) = C_i(Ω, Δ)^{-1} (2π)^{-J/2} |Ω|^{-1/2} exp[-1/2(θ_i - Δz_i)'Ω^{-1}(θ_i - Δz_i)]
```

where `C_i(Ω, Δ) = P(u_i ≤ -Δz_i | u_i ~ N_J(0, Ω))` is the integrating constant (a J-dimensional multivariate normal integral).

**This normalizing constant is data-dependent (depends on z_i) and is NOT included in the current Stan code!**

**Correction needed:** The log probability must include `-log(C_i(Ω, Δ))` for each observation. This is a multivariate normal CDF that needs to be computed.

---

### 2. CRITICAL: Technical Inefficiency u₀ Distribution is INCORRECTLY Specified

**Paper specification (Assumption 2, page 5):**

```
u₀_i ~ iid N(0, σ²_u), u₀_i ≤ 0
```

This is a **half-normal distribution** truncated at 0 from above (negative half-normal).

**Current implementation in `gsf_model.stan` (line 92):**

```stan
u0 ~ normal(0, sigma_u);
```

**Problem:** Again, this uses an ordinary normal, not a truncated normal. While the parameter is constrained by `vector<upper=0>[N] u0;`, the **normalizing constant** `(π σ²_u / 2)^{-1/2}` (see equation 8 in the paper) is missing.

**Correction needed:** For a half-normal (normal truncated at 0):

```stan
u0 ~ normal(0, sigma_u) T[, 0];  // Truncated normal with upper bound at 0
```

Or equivalently, add the log of the normalizing constant:

```stan
target += normal_lpdf(u0 | 0, sigma_u) - N * normal_lccdf(0 | 0, sigma_u);
```

---

### 3. ISSUE: Elasticity Calculation in Output Loss Uses Wrong Inputs

**Paper specification (Equation 5, page 5):**

The composed error term is:

```
v_i = [θ'_i q_i + (1/2)θ'_i B θ_i] + u₀_i + ε₀_i
```

where `q_i = α + Bx_i` (the elasticity at the **observed** inputs, not effective inputs).

**Current implementation in `gsf_model.stan` (lines 143-144, 150-152):**

```stan
// Output elasticities: q = alpha + B*x_eff
q_n = alpha + B * x_eff;
...
// Output loss from slacks: theta'q + 0.5*theta'B*theta
output_loss_from_slacks[n] = -(dot_product(theta_n, q_n) + 0.5 * quad_form(B, theta_n)) * 100;
```

**Problem:** The paper defines `q_i = α + Bx_i` (using **observed** inputs x_i), but the code uses `q_n = alpha + B * x_eff` (using **effective** inputs x_eff = x_i + θ_i).

Looking more carefully at equation (5), we have:

```
[θ'_i α + (1/2)θ'_i Bθ_i) + θ'_i Bx_i]
```

which can be rewritten as:

```
θ'_i(α + Bx_i) + (1/2)θ'_i Bθ_i = θ'_i q_i + (1/2)θ'_i Bθ_i
```

where `q_i = α + Bx_i`. **The elasticity in this context should use x_i, NOT x_eff!**

**Correction needed:** Change:

```stan
q_n = alpha + B * x_n;  // Use observed inputs, not effective inputs
```

---

### 4. Elasticity Reporting (Minor Clarification)

**Paper specification (Part 2 of Appendix, and equation in Section 5 page 8):**

```
Elas_i = ∂E(y_i|z_i,x_i)/∂x_i = α + B(x_i + θ̄_i)
```

where `θ̄_i = E(θ_i|y,X,Z)` is the **posterior mean** of θ.

**Current implementation in `gsf_model.stan` (lines 143-145):**

```stan
// Output elasticities: q = alpha + B*x_eff
q_n = alpha + B * x_eff;
```

Here `x_eff = x_n + theta_n` where `theta_n` is the MCMC draw.

**This is correct** for generating posterior samples of elasticities that are then averaged in R. The elasticities reported in Table 2 use effective inputs (x + θ), which the current code does correctly.

**Note:** There are TWO different elasticity concepts in the paper:
1. For **output loss calculation** (equation 5): use `q_i = α + Bx_i` (observed inputs)
2. For **elasticity reporting** (Table 2): use `α + B(x_i + θ_i)` (effective inputs)

The current code conflates these two uses.

---

### 5. Technical Inefficiency Percentage Calculation (Acceptable)

**Paper Table 2:** Reports mean technical inefficiency as 2.43% (with endogeneity correction).

**Current implementation in `gsf_model.stan` (line 136):**

```stan
technical_inefficiency[n] = (exp(-u0[n]) - 1) * 100;
```

Since u₀ ≤ 0, we have -u₀ ≥ 0, so exp(-u₀) ≥ 1, and (exp(-u₀) - 1) ≥ 0. This gives the percentage output loss.

The paper's Table 2 reports "u₀ (technical inefficiency)" with mean 0.024 (2.4%), suggesting they may report -u₀ directly. However, for small values, `exp(-u₀) - 1 ≈ -u₀`, so the difference is negligible.

**This is acceptable.**

---

### 6. Input Slack Percentage Calculation (Acceptable)

**Paper interpretation:** θ_j ≤ 0 means `-θ_j × 100` gives the percentage over-use of input j.

**Current implementation (line 140):**

```stan
input_slacks[n, j] = (exp(-theta[n, j]) - 1) * 100;
```

This uses the exact formula for percentage. For small θ, this approximates `-θ × 100`.

**This is acceptable.**

---

### 7. Issue in Alternative Stan File (`src/stan_files/gsf.stan`)

The alternative Stan file has similar issues plus additional problems:

- Uses `cov_matrix[J] Omega` directly which is less efficient than Cholesky parameterization
- Uses `inv_gamma(1, 1)` priors for variances instead of proper half-Cauchy or exponential
- No truncation correction for either θ or u₀
- Same conflation of elasticity concepts

---

## Simulation Data Generation Review

The `simulate_gsf_data.R` and `simulate_uk_manufacturing.R` functions **correctly** generate data according to the model:

- θ is truncated at 0 from above (line 162: `pmin(theta_i, -0.001)`)
- u₀ is generated as negative half-normal (line 166: `-abs(rnorm(N, 0, sigma_u))`)
- The translog production function is correctly specified

This means the **data generation is correct**, but the **estimation** (Stan model) does not match the data generating process exactly due to missing truncation corrections.

---

## Summary of Required Corrections

### Critical Fixes:

#### 1. Fix truncated normal for u₀

```stan
// In model block, replace:
u0 ~ normal(0, sigma_u);

// With:
for (n in 1:N) {
    u0[n] ~ normal(0, sigma_u) T[, 0];
}

// Or add the normalizing constant manually:
target += normal_lpdf(u0 | 0, sigma_u);
target += -N * log(2);  // Correction for half-normal (integral from -∞ to 0 is 0.5)
```

#### 2. Fix truncated multivariate normal for θ

This is more complex. The paper's MCMC approach (Appendix B) uses a Gibbs sampler that handles truncated MVN directly using specialized algorithms (references Geweke 1991).

In Stan, you have several options:

**Option A: Add the log-normalizing constant explicitly**

```stan
// This requires computing the MVN CDF, which Stan doesn't have built-in
// For J=2 (bivariate case), you can use:
functions {
    real binormal_cdf(real x1, real x2, real rho) {
        // Implementation of bivariate normal CDF
        // See Drezner (1992) or use numerical integration
    }
}

model {
    for (n in 1:N) {
        vector[J] mu_theta = Delta * Z_slack[n]';
        target += multi_normal_cholesky_lpdf(theta[n]' | mu_theta, L_Sigma);

        // Subtract log of truncation probability
        // For J=2: use bivariate normal CDF
        real C_i = binormal_cdf(-mu_theta[1]/sigma_slack[1],
                                 -mu_theta[2]/sigma_slack[2],
                                 L_Omega[2,1]);
        target += -log(C_i);
    }
}
```

**Option B: Use a different parameterization**

Transform θ to unconstrained space using a bijective mapping (e.g., log transform for negative values), then adjust the Jacobian.

**Option C: Accept the approximation**

The current approach with parameter constraints provides an approximation. Document this limitation and note that for observations where the mean Δz_i is far from 0, the bias may be larger.

#### 3. Fix the Output Loss Calculation

In `gsf_model.stan`, in the `generated quantities` block:

```stan
// For output loss calculation, use observed inputs:
vector[J] q_for_loss = alpha + B * x_n;  // NOT x_eff
output_loss_from_slacks[n] = -(dot_product(theta_n, q_for_loss) +
                                0.5 * quad_form(B, theta_n)) * 100;

// For elasticity reporting, continue using effective inputs:
vector[J] q_n = alpha + B * x_eff;
elasticities[n] = q_n';
```

---

## Impact Assessment

### Without the truncation corrections:

The model will sample from the **wrong posterior distribution**. Parameter estimates will be biased because the likelihood is not correctly specified. The severity depends on:

- How far Δz_i is from 0 (if means are strongly negative, truncation matters less)
- The variance of the slacks (smaller variance = truncation matters more)

### Without the elasticity/output loss fix:

Output loss calculations will be incorrect. This affects post-estimation interpretation but not parameter estimation.

---

## Technical Note on Identification

The paper emphasizes that identification comes from:

1. θ appearing both linearly (θ'α) and quadratically ((1/2)θ'Bθ)
2. θ interacting with x (θ'Bx)
3. Distributional assumptions on θ, u₀, and ε

The **translog specification is essential** for identification - a Cobb-Douglas production function would not identify the input-specific slacks without additional structure (like input-specific z variables as in KOS).

---

## Recommendations

### For an Exact Implementation:

1. **Implement proper truncation for u₀** using Stan's `T[, 0]` syntax
2. **Implement proper truncation for θ** by adding the multivariate normal CDF normalizing constant
3. **Separate the two elasticity calculations** (for output loss vs. for reporting)

### For a Pragmatic Implementation:

1. At minimum, fix the u₀ truncation (easy in Stan)
2. Document the θ truncation approximation as a limitation
3. Fix the output loss calculation to use observed inputs

### For Validation:

1. Run simulations comparing true parameters to estimated parameters
2. The `gsf_recovery` dataset or similar simulation studies can verify parameter recovery
3. Compare efficiency estimates to true values in simulated data

---

## Files Reviewed

| File | Purpose | Issues Found |
|------|---------|--------------|
| `inst/stan/gsf_model.stan` | Main Stan model | Critical: truncation, output loss |
| `src/stan_files/gsf.stan` | Alternative Stan model | Same issues + inefficient priors |
| `R/gsf_fit.R` | Model fitting | OK |
| `R/efficiency.R` | Efficiency extraction | OK (uses Stan output) |
| `R/simulate_data.R` | Data simulation | Correct |
| `R/summary.R` | Summary methods | OK |

---

## References

- Kumbhakar, S.C., Tsionas, M.G. (2021). "Dissections of input and output efficiency: A generalized stochastic frontier model." *International Journal of Production Economics*, 232, 107940.
- Geweke, J. (1991). "Efficient simulation from the multivariate normal and student-t distributions subject to linear constraints."
- Drezner, Z. (1992). "Computation of the multivariate normal integral." *ACM Transactions on Mathematical Software*, 18, 470-480.
