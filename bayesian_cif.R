library(gtools)  # For rdirichlet()

#' Compute Posterior Samples of the Bayesian Cumulative Incidence Function (CIF)
#'
#' @param time Vector of time-to-event values.
#' @param event1 Binary vector for event type 1 (e.g., CV death).
#' @param event2 Binary vector for event type 2 (e.g., non-CV death).
#' @param alpha_prior Dirichlet prior vector of length 3.
#' @param posterior_sample_size Number of posterior samples to draw.
#' @param event_index Integer: 1 for event1, 2 for event2.
#' @param interval Time grid spacing (default = 1).
#'
#' @return A data.frame of posterior CIF samples:
#'         rows = time points, columns = posterior samples.
bayesian_cif <- function(time,
                         event1,
                         event2,
                         alpha_prior = c(0.0001, 0.0001, 0.0001),
                         posterior_sample_size = 100000,
                         event_index = 1,
                         interval = 1) {

  # Input validation
  if (!event_index %in% c(1, 2)) {
    stop("event_index must be either 1 (for event1) or 2 (for event2).")
  }

  # Create discretized time grid
  discretized_time <- seq(0, max(time), by = interval)

  # Initialize variables
  n <- length(time)
  current_risk_set <- n
  product <- matrix(1, nrow = posterior_sample_size, ncol = 1)
  cumulative_incidence <- matrix(0, nrow = length(discretized_time), ncol = posterior_sample_size)

  # Loop over time points (except last)
  for (t in 1:(length(discretized_time) - 1)) {

    # Event counts at current time
    count_e1 <- sum(event1[time == discretized_time[t]])
    count_e2 <- sum(event2[time == discretized_time[t]])

    # Posterior sampling
    observed_counts <- c(count_e1, count_e2, current_risk_set)
    posterior_samples <- rdirichlet(posterior_sample_size, observed_counts + alpha_prior)

    # Update CIF for specified event index
    cumulative_incidence[t + 1, ] <- cumulative_incidence[t, ] + posterior_samples[, event_index] * product

    # Update survival probability
    product <- product * (1 - posterior_samples[, 1] - posterior_samples[, 2])

    # Update risk set
    current_risk_set <- sum(time >= discretized_time[t])
  }

  # Return posterior samples as data.frame
  cif_df <- as.data.frame(cumulative_incidence, row.names = NULL)
  return(cif_df)
}
