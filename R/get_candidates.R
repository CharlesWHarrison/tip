#' @title Compute the set of similar subjects
#' @description Find the <.num_candidates> subjects that are most similar to subject .i.
#' @param .i Positive integer: a positive integer corresponding to the subject index (e.g., a row index in a data frame for vector-variate data).
#' @param .similarity_matrix Matrix: an \code{n} x \code{n} matrix of similarity values.
#' @param .num_candidates Positive integer: the number of similar subjects to be extracted.
#' @returns Vector of positive integers: a vector of the <\code{.num_candidates}> subject indices that are most similar to the subject with index <\code{.i}>.
#' @noRd
get_candidates <- function(.i, .similarity_matrix, .num_candidates){
  # --- A function to return the .num_candidates indices corresponding to
  # the subjects that are most similar to subject .i ---
  # Note: start at 2 since 1 is always the candidate itself
  return(order(.similarity_matrix[.i,], decreasing = TRUE)[2:(.num_candidates + 1)])
}
