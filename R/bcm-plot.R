#' @title Generate plots from a Bayesian Clustering Model (bcm) object
#' @rdname plot.bcm
#' @aliases plot
#' @param x bcm object: a Bayesian Clustering Model (bcm) object fit to a dataset
#' @param y Not used.
#' @param ... Not used.
#' @example example/bcm-plot_examples.R
#' @exportMethod plot
setMethod("plot", signature(x="bcm",y="missing"), function(x,y,...){
  return(list(trace_plot_posterior_number_of_clusters = ggplot_number_of_clusters_trace(.posterior_number_of_clusters = x@posterior_number_of_clusters),
              histogram_posterior_number_of_clusters = ggplot_number_of_clusters_hist(.posterior_number_of_clusters = x@posterior_number_of_clusters)))
}
)
