#' @title Compute a proximity matrix
#' @description A function to convert a vector of posterior cluster assignments into
#' an n x n matrix B where Bij = 1 if vector[i] == vector[j] and 0 otherwise. That is, Bij = 1 if
#' the (i)th subject is in the same cluster as the (j)th subject and Bij = 0 if the (i)th subject
#' and the (j)th subject are not in the same cluster.
#' @param .assignments Vector of positive integers: the (i)th element denotes the (i)th subject's cluster assignment after posterior sampling.
#' @returns Matrix: an \code{n} x \code{n} matrix B where Bij = 1 if vector[i] == vector[j] and 0 otherwise (i.e., Bij = 1 if the (i)th subject and the (j)th subject are in the same cluster and Bij = 0 if the (i)th subject and the (j)th subject are not in the same cluster).
#' @noRd
get_proximity_matrix <- function(.assignments){
  # --- A function to construct a proximity matrix based on
  # a vector of .assignments ---
  # matrix_{i,j} = 1 if subject_i and subject_j belong to the same cluster
  # matrix_{i,j} = 0 if subject_i and subject_j do not belong to the same cluster
  return(outer(.assignments, .assignments, function(x, y) as.integer(x==y)))
}
