# Bayesian Cumulative Incidence Function (CIF) Tutorial

This repository provides an R function to compute posterior samples of the **Bayesian cumulative incidence function (CIF)** for competing risks using a **Dirichlet–multinomial** approach.  

The method is fully nonparametric and does not rely on proportional hazards or cause-specific hazard models. Instead, it uses a discrete-time Bayesian formulation with a Dirichlet prior on the daily hazard probabilities, producing posterior draws of the CIF that can be directly summarized and visualized.  

This implementation accompanies the manuscript:  
> **Daniel Paydarfar, Brian L. Claggett, Lu Tian, Scott D. Solomon, and Lee-Jen Wei. “Bayesian Survival Analysis for Comparative Cardiovascular Clinical Trials with Competing Risks: Application to the PARAGON-HF Study”** 

--------------------------------------------------------------------------------
## 1. Overview

### 1.1 What Is `bayesian_cif`?

The `bayesian_cif` function generates **posterior draws** of the CIF for a specified event type in the presence of competing risks. Each time point is treated as a multinomial trial with three possible outcomes:  

1. Event type 1 (e.g., CV death)  
2. Event type 2 (e.g., non-CV death)  
3. No event (survival to the next interval)  

A Dirichlet prior ensures conjugacy, so the posterior is also Dirichlet and can be sampled efficiently.

- **Inputs**:  
  - `time`: numeric vector of observed times (event or censoring).  
  - `event1`: binary indicator vector for event type 1.  
  - `event2`: binary indicator vector for event type 2.  
  - `event_index`: choose which event’s CIF to compute (`1` for event1, `2` for event2).  
  - `alpha_prior`: length-3 vector of Dirichlet prior parameters (default: `c(0.0001, 0.0001, 1)`).  
  - `posterior_sample_size`: number of posterior samples (default: 100000).  
  - `interval`: discretization interval for time (default: 1).  

- **Output**:  
  A matrix whose rows correspond to discrete time points (from 0 to max observed time, spaced by `interval`), and whose columns correspond to posterior samples of the CIF trajectory for the chosen event.  

### 1.2 Why Use This Approach?

- **Model-free**: No assumptions of proportional hazards.  
- **Fully Bayesian**: Uncertainty is captured naturally via posterior samples, yielding interpretations like "there’s a 95% chance the true survival or RMST lies in this range”.  
- **Flexible**: Handles arbitrary numbers of competing risks with customizable priors.  
- **Efficient**: Conjugacy of the Dirichlet–multinomial ensures fast computation.  

--------------------------------------------------------------------------------
## 2. Installing / Loading the Function

Load the function directly into R via the following command in R:

```r
# Source directly from GitHub
source("https://raw.githubusercontent.com/danali-bst/bayesian-cif/main/bayesian_cif.R")
```

Now `bayesian_cif` is available in your R session.

--------------------------------------------------------------------------------
## 3. Example: Two Competing Risks

Suppose we have a dataset with two event types: cardiovascular (CV) death and non-CV death.  

```r
# Example data: 12 patients
df <- data.frame(
  time     = c(1, 2, 2, 3, 5, 5, 7, 8, 10, 10, 12, 15),
  event_cv = c(1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0),  # Event type 1
  event_ncv= c(0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1)   # Event type 2
)

# Run Bayesian CIF for event type 1 (CV death)
set.seed(123)
cif_samples <- bayesian_cif(
  time = df$time,
  event1 = df$event_cv,
  event2 = df$event_ncv,
  event_index = 1,             # compute CIF for CV death
  alpha_prior = c(0.0001, 0.0001, 1),
  posterior_sample_size = 5000,
  interval = 1
)

dim(cif_samples)
# rows = number of discrete time points, cols = posterior samples
```

### 3.1 Summarizing the Posterior (e.g. Posterior Mean)

```r
# Posterior mean CIF
mean_cif <- apply(cif_samples, 1, mean)

# 95% credible interval at each time
lower_ci <- apply(cif_samples, 1, quantile, 0.025)
upper_ci <- apply(cif_samples, 1, quantile, 0.975)

# Plot
time_grid <- seq(0, max(df$time), by = 1)
plot(time_grid, mean_cif, type = "s",
     ylim = c(0, 1),
     xlab = "Time",
     ylab = "Cumulative Incidence (Event 1)",
     main = "Bayesian CIF with 95% Credible Interval")
lines(time_grid, lower_ci, type = "s", lty = 2, col = "red")
lines(time_grid, upper_ci, type = "s", lty = 2, col = "red")
```

Interpretation: The dashed bands represent posterior uncertainty in the CIF trajectory.

--------------------------------------------------------------------------------
## 4. Additional Summaries

### Cause-Specific Restricted Mean Time Lost (RMTL)

The area under the CIF curve up to time `tau` is the **restricted mean time lost** due to that cause (complementary to RMST).  

```r
# Compute RMTL up to tau
tau = 10
rmtl_draws <- colSums(cif_samples[1:tau, ])
rmtl_mean  <- mean(rmtl_draws)
rmtl_ci    <- quantile(rmtl_draws, c(0.025, 0.975))

cat("Posterior mean RMTL up to day", tau, "=", rmtl_mean, "\n")
cat("95% Credible Interval =", rmtl_ci, "\n")
```

This corresponds to the “AUC of the CIF” discussed in the manuscript.

--------------------------------------------------------------------------------
## 5. Extensions

- Can be generalized to >2 competing risks by extending the Dirichlet dimension.  
- Priors can be adjusted to reflect prior knowledge about relative event frequencies.  
- Posterior samples can be combined across strata (similar to the [Bayesian survival tutorial](https://github.com/danali-bst/bayesian-survival-tutorial)).  

--------------------------------------------------------------------------------
## 6. Closing Notes

This Bayesian CIF estimator provides a simple and flexible alternative to cause-specific or subdistribution hazards. It is especially useful for visualizing uncertainty and performing inference without relying on hazard models.  

For questions or contributions, please open an issue or pull request in this repository, or contact Daniel Paydarfar at danielpaydarfar@fas.harvard.edu.

--------------------------------------------------------------------------------
## 7. Appendix (mathematical description of method)

Let us consider a scenario involving two types of events: first heart failure hospitalization or CV death and non-CV death as in the PARAGON-HF trial. Extending this case to more than two events of interest is straightforward.

Now, in practice, "day" is typically used as the smallest time unit for recording survival events. Accordingly, we consider analysis based on the data collected at discrete time points:

$$
t = 1, \ldots, M.
$$

For example, in the PARAGON-HF data, \( M = 1500 \). On day \( t \), let \( r_1(t) \) represent the cause-specific hazard of first heart failure hospitalization or CV death, which is the conditional probability that a patient experiences the composite outcome on day \( t \), given survival up to that day. Similarly, let \( r_2(t) \) be the conditional probability of non-CV death. The overall survival function at day \( t \) can be expressed as

$$
S(t) = Pr\ (survival\ time\ >\ t) = \prod_{s=1}^{t}\left(1-r_1(s)-r_2(s)\right).
$$

At each time \( t \), a patient's outcome follows a three-category multinomial distribution with two free parameters, \( r_1(t) \) and \( r_2(t) \). We consider a conjugate prior for this multinomial variable: a Dirichlet distribution with parameters \( (\alpha_1,\alpha_2,\alpha_3) \). For the entire \( M \) days, the prior is the product of \( M \) individual Dirichlet priors, which is proportional to

$$
\prod_{t=1}^{M} r_1(t)^{\alpha_1} r_2(t)^{\alpha_2} \left(1-r_1(t)-r_2(t)\right)^{\alpha_3}
$$

To construct the likelihood function for the parameters \( \{r_1(t),\ r_2(t)\colon\ t=1,\ \ldots,\ M\} \), the contribution from an individual patient observed at time \( t \) is as follows:

- If the patient experiences the composite outcome (category 1), the contribution is \( r_1(t)\times S(t-1) \)
- If the patient dies due to non-CV causes (category 2), the contribution is \( r_2(t)\times S(t-1) \)
- If the patient is alive after \( t \) (category 3), the contribution is \( (1-r_1(t)-r_2(t))\times S(t-1) \)

The full likelihood \( L \) is proportional to the product of these individual contributions across \( M \) days:

$$
\prod_{t=1}^{M} r_1(t)^{d_1(t)} r_2(t)^{d_2(t)} \left(1-r_1(t)-r_2(t)\right)^{n(t)-d_1(t)-d_2(t)},
$$

where \( d_1(t) \) and \( d_2(t) \) are the observed numbers of patients in categories 1 and 2, respectively, and \( n(t) \) is the number of patients who survived past \( t-1 \). The likelihood function \( L \) can be interpreted similarly to treating \( M \) independent multinomial experiments, a result justified through survival analysis martingale theory [14].

At time \( t \), the posterior distribution is also Dirichlet with parameters \( (\alpha_1+d_1(t),\ \alpha_2+d_2(t),\ \alpha_3+n(t)-d_1(t)-d_2(t)) \). The complete posterior distribution is the product of these \( M \) individual posterior distributions. One can then generate a large number of realizations of this Dirichlet distribution to approximate the posterior distribution. The cumulative incidence function for first heart failure hospitalization or CV death is given by:

$$
F_1(t)=\sum_{s=1}^{t} r_1(s)\prod_{u=1}^{s-1}\left(1-r_1(u)-r_2(u)\right).
$$

The posterior distribution of CIF can be obtained through the previous realizations of \( \{r_1(t),\ r_2(t)\colon\ t=1,\ \ldots,\ M\} \). Lastly, the area under the \( F_1(t) \) can be approximated and used to compare two treatment groups by examining the difference or ratio of their AUCs.

When there is limited information about the daily hazards, a common practice is to choose a prior which would have a small impact on the final analysis. That is, we may choose very small parameter values for \( (\alpha_1,\alpha_2,\alpha_3) \).



