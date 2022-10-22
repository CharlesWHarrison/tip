#' @title Recode a vector of positive integer values to start at 1
#' @description A function to "recode" a vector. For example, 2, 3, 5, 6, 10, 2, 2, 2, 5 needs
#' to be recoded to 1, 2, 3, 4, 5, 1, 1, 1, 3. This function is used to ensure that the posterior
#' cluster assignments start at 1 (otherwise an error occurs).
#' @param .posterior_assignments Vector of positive integers; each positive integer corresponds to a posterior cluster assignment.
#' @noRd
recode <- function(.posterior_assignments){
  # --- A function to recode the current cluster assignments so that each
  # cluster assignment is in the set {1, 2, 3, ..., K} and there are no gaps.
  # Example: 2, 3, 5, 6, 10, 2, 2, 2, 5 needs to be recoded to 1, 2, 3, 4, 5, 1, 1, 1, 3
  # Example: 1, 2, 3, 4, 5, 2, 2, 2, 5 is recoded to itself since it is already correct ---

  # Compute the current sorted unique cluster assignment values
  .sorted_unique_cluster_values <- sort(unique(.posterior_assignments), decreasing = FALSE)

  # For each current cluster assignment value
  for(j in 1:length(.posterior_assignments)){
    # Replace the current cluster assignment value with the index corresponding to the
    # sorted unique cluster value that the current cluster assignment is equal to
    .posterior_assignments[j] <- which(.sorted_unique_cluster_values == .posterior_assignments[j])
  }
  return(.posterior_assignments)
  # recode(.posterior_assignments = c(2, 3, 5, 6, 10, 2, 2, 2, 5))
}
