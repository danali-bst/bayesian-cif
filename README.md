# Bayesian Cumulative Incidence Function estimator in R

An R function to compute posterior samples of the Bayesian cumulative incidence function (CIF) for competing risks using a Dirichlet-multinomial model.

# Description

This function returns posterior sample trajectories for the CIF of a specified event type in the presence of competing risks. It uses a nonparametric Bayesian approach with a Dirichlet prior and allows fine control over the number of posterior samples and time discretization.

# Installation

No package installation is needed — you can source the function directly into R using:

source(“https://raw.githubusercontent.com/YOUR_USERNAME/bayesian-cif/main/bayesian_cif.R”)

# Usage

cif_samples <- bayesian_cif(
time = time_vector,
event1 = event_cv,      # binary vector for event type 1
event2 = event_ncv,     # binary vector for event type 2
event_index = 1,        # choose 1 for event1 or 2 for event2
alpha_prior = c(0.0001, 0.0001, 0.0001),  # optional
posterior_sample_size = 10000,           # optional
interval = 1                             # optional
)

Arguments
	•	time: Numeric vector of time-to-event or censoring values.
	•	event1: Binary vector for event type 1 (e.g., CV death).
	•	event2: Binary vector for event type 2 (e.g., non-CV death).
	•	alpha_prior: Vector of length 3 with Dirichlet prior parameters (default: c(0.0001, 0.0001, 0.0001)).
	•	posterior_sample_size: Number of posterior samples to generate (default: 100000).
	•	event_index: Integer: 1 for event1, 2 for event2.
	•	interval: Time discretization interval (default: 1).

Output

A data.frame with:
	•	Rows = time points (from 0 to max observed time, spaced by interval)
	•	Columns = posterior samples of the CIF for the chosen event

Each column represents a full posterior trajectory for the CIF. You can compute point estimates, credible intervals, or plot directly from this matrix.

# Example

cif_samples <- bayesian_cif(
time = df$time,
event1 = df$event_cv,
event2 = df$event_ncv,
event_index = 1,
posterior_sample_size = 5000,
interval = 10
)

Summarize or visualize

mean_cif <- apply(cif_samples, 1, mean)
