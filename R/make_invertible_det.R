#' @title Automatic matrix inversion checking via a determinant criterion
#' @param .matrix Square matrix: the matrix that may or may not be invertible.
#' @param .step_size Positive numeric: if necessary, a small value that is ITERATIVELY added to each diagonal element of a matrix until the matrix is invertible.
#' @returns Matrix: a matrix that can be inverted.
#' @noRd
make_invertible_det <- function(.matrix, .step_size = 0.01){
  # --- A function used to make a matrix invertible by
  # adding a small number to the matrix's diagonal.
  # If the matrix is invertible, then return the matrix. ---
  if(det(.matrix)< 0.001){
    .bias_size <- .step_size
    num_rows_cols <- dim(.matrix)[1]
    .temp <- .matrix + diag(num_rows_cols)*.bias_size
    if(det(.temp) > 0.001){.matrix = .temp}
    while(det(.temp) <= 0.001){
      .temp = .matrix + diag(num_rows_cols)*.bias_size
      .bias_size = .bias_size + .step_size
    }
    .matrix <- .temp
  }
  return(.matrix)
}
