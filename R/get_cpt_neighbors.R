#' @title Estimate The Number of Similar Subjects
#' @description Estimate the number of similar subjects using univariate multiple change point detection (i.e. binary segmentation).
#' @param .distance_matrix Matrix. A symmetric n x n matrix of distance values.
#' @returns Vector of positive integers. A vector of positive integers where the (i)th integer corresponds to the number of subjects (observations) that are similar to the (i)th subject (observation).
#' @examples
#' X <- matrix(rnorm(10*10),nrow=20,ncol=5)
#' get_cpt_neighbors(.distance_matrix = data.matrix(dist(X)))
#' @export
get_cpt_neighbors <- function(.distance_matrix){
  # --- A function used to obtain the nearest
  # neighbors for each subjects based on their
  # mutual distances ---

  # A vector to hold the number of subjects
  # most similar to subject i
  .num_neighbors <- vector()

  # Compute the number of subjects
  .n <- dim(.distance_matrix)[1]
  # Apply change-point detection to the set of
  # sorted distances corresponding to each subject i
  for(i in 1:.n){
    .num_neighbors[i] <- changepoint::cpt.mean(data = sort((.distance_matrix[i,])[-i],
                                                           decreasing = FALSE),
                                               method = "BinSeg", Q = floor(.n/2))@cpts[1]
  }
  return(.num_neighbors)
  # get_cpt_neighbors(.distance_matrix = matrix(abs(rnorm(100)),nrow=10,ncol=10))
}