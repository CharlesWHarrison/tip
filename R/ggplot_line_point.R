#' @title Plot connected points using ggplot2
#' @description A function to that produces a ggplot2 plot of .y versus .x
#' where points are added via geom_point() and the points are connected via geom_line().
#' @param .x The variable on the horizontal axis.
#' @param .y The variable on the vertical axis.
#' @param .ylab Character; the label on the vertical axis.
#' @param .xlab Character; the label on the horizontal axis.
#' @returns ggplot2 plot; a ggplot2 plot of .y versus .x with a label .xlab on the horizontal axis and label .ylab on the vertical axis.
#' @importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab
#' @importFrom rlang .data
#' @examples ggplot_line_point(.x = rnorm(10), .y = rnorm(10), .xlab = "x", .ylab = "y")
#' @export
ggplot_line_point <- function(.x, .y, .xlab = "", .ylab = ""){
  # --- A function to plot a line and the corresponding points ---
  plot_df <- data.frame(x = .x, y = .y)
  p <- ggplot(data = plot_df, aes(x = .data$x, y = .data$y)) + geom_line()
  p <- p + geom_point() + xlab(.xlab) + ylab(.ylab)
  return(p)
}
