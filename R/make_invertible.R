#' @title Automatic matrix inversion checking via an eigenvalue criterion
#' @param .matrix Square matrix; the square matrix that may or may not be invertible.
#' @param .step_size Positive numeric; if necessary, a small value that is ITERATIVELY added to each diagonal element of a matrix until the matrix is invertible.
#' @noRd
make_invertible <- function(.matrix, .step_size = 0.01){
  # --- A function used to make a matrix invertible by
  # adding a small number to the matrix's diagonal.
  # If the matrix is invertible, then return the matrix. ---
  .bias_size <- .step_size
  # Compute the number of columns
  .dimension <- dim(.matrix)[2]
  if(any(eigen(.matrix)$values <= 0.0)){
    .temp <- .matrix + diag(.dimension)*.bias_size
    if(all(eigen(.matrix)$values > .step_size)){return(.temp)}
    while(any(eigen(.temp)$values < 0.0)){
      .temp = .matrix + diag(.dimension)*.bias_size
      .bias_size = .bias_size + .step_size
    }
    .matrix <- .temp
  }
  return(.matrix)
  # make_invertible(.matrix = matrix(data = 1, nrow = 3, ncol = 3),
  #                 .step_size = 0.01)
}
