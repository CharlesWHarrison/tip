#' @title Compute the prior probability using the Table Invitation Prior (TIP)
#' @description Compute the prior probability that a subject belongs to a cluster.
#' @param .i Positive integer: the subject index (i.e. row index in a data frame for vector-variate data or the (.i)th matrix or (.i)th tensor).
#' @param .similarity_matrix Matrix: a matrix of pairwise subject similarity values.
#' @param .current_assignments Vector of integers: each integer is the posterior cluster assignment after the invitation step.
#' @param .num_clusters Positive integer: the number of clusters after the invitation step.
#' @returns Vector: a vector of prior probabilities that the (.i)th subject belongs to each cluster.
#' @noRd
prob_tip_i <- function(.i, .similarity_matrix, .current_assignments, .num_clusters){
  # --- A function to compute the conditional posterior probability of a
  # subject joining a cluster (table) ---

  # Define a vector to hold the probability for subject i to be in cluster k = 1, 2, ..., <.num_clusters>
  .prob_k_vector <- vector()

  # Compute sum lambda(i,j) for j s.t. cluster[j] == k
  for(.k in 1:.num_clusters){
    .prob_k_vector[.k] <- sum(.similarity_matrix[.i, which(.current_assignments == .k)])
  }

  # Return the vector of probabilities
  return(.prob_k_vector)
}
