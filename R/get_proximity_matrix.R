#' @title Compute a proximity matrix
#' @description A function to convert a vector of posterior cluster assignments into
#' an n x n matrix B where Bij = 1 if vector[i] == vector[j] and 0 otherwise.
#' @param .assignments Vector of positive integers; the (i)th element denotes the (i)th subjects cluster assignment after posterior sampling.
#' @noRd
get_proximity_matrix <- function(.assignments){
  # --- A function to construct a proximity matrix based on
  # a vector of .assignments ---
  # matrix_{i,j} = 1 if subject_i and subject_j belong to the same cluster
  # matrix_{i,j} = 0 if subject_i and subject_j do not belong to the same cluster
  return(outer(.assignments, .assignments, function(x, y) as.integer(x==y)))
}
