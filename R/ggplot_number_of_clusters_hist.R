#' @title Plot the posterior distribution of the number of clusters.
#' @description A function that produces a ggplot bar chart (i.e., geom_bar) that corresponds
#' to the posterior number of clusters. The vertical axis is normalized so that it displays
#' the posterior probability.
#' @param .posterior_number_of_clusters Vector of positive integers: each integers corresponds to the
#'  number of clusters after posterior sampling for each iteration in the Gibbs sampler.
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab scale_x_continuous
#' @importFrom rlang .data
#' @returns ggplot2 geom_bar plot: a plot of the distribution of the posterior number of clusters
#' computed after each sampling iteration in the Gibbs sampler.
#' @example example/ggplot_number_of_clusters_hist_examples.R
#' @export
ggplot_number_of_clusters_hist <- function(.posterior_number_of_clusters){
  .table = table(.posterior_number_of_clusters)
  # --- A function to construct the posterior distribution of the number of clusters ---
  .df_tip <- data.frame(y = as.numeric(.table)/sum(.table),
                        x = as.numeric(names(.table)))
  plot <- ggplot(data = .df_tip) + geom_bar(aes(x = .data$x, y = .data$y), stat = "identity")
  plot <- plot + xlab("Number of Clusters")
  plot <- plot + ylab("Posterior Probability")
  plot <- plot + scale_x_continuous(breaks = 1:max(.df_tip$x))
  return(plot)
}
