#' @title Compute the set of similar subjects
#' @description Find the <.num_candidates> subjects that are most similar to subject .i.
#' @param .i Positive integer; the subject index (i.e. row index in a matrix for vector-variate data).
#' @param .similarity_matrix Matrix; an n x n matrix of similarity values.
#' @param .num_candidates Positive integer; the number of similar subjects to be extracted.
#' @noRd
get_candidates <- function(.i, .similarity_matrix, .num_candidates){
  # --- A function to return the .num_candidates indices corresponding to
  # the subjects that are most similar to subject .i ---
  # Note: start at 2 since 1 is always the candidate itself
  return(order(.similarity_matrix[.i,], decreasing = TRUE)[2:(.num_candidates + 1)])
}
