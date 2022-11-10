#' @title Plot the trace plot of the posterior number of clusters
#' @description A function that produces a ggplot2 trace plot (i.e., geom_line)
#' with respect to the posterior number of clusters.
#' @param .posterior_number_of_clusters Vector of positive integers: each (s)th element denotes the number of clusters after posterior sampling
#' for each iteration s = 1, 2, ..., \code{samples} + \code{burn} in the Gibbs sampler.
#' @returns ggplot2 geom_line plot: a plot of the posterior number of clusters in each Gibbs sampling iteration versus the Gibbs sampling iteration number.
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab
#' @importFrom rlang .data
#' @example man/example/ggplot_number_of_clusters_trace_examples.R
#' @export
ggplot_number_of_clusters_trace <- function(.posterior_number_of_clusters){
  # --- A function to construct the posterior distribution of the number of clusters ---
  .df_tip <- data.frame(iteration = 1:length(.posterior_number_of_clusters),
                        num_clusters = .posterior_number_of_clusters)
  plot <- ggplot(data = .df_tip, aes(y = .data$num_clusters, x = .data$iteration)) + geom_line()
  plot <- plot + xlab("Gibbs Sampling Iteration")
  plot <- plot + ylab("Posterior Number of Clusters")
  return(plot)
}
