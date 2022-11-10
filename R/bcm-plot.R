#' @title Generate plots from a Bayesian Clustering Model (bcm) object
#' @rdname plot.bcm
#' @aliases plot
#' @param x bcm object: a Bayesian Clustering Model (bcm) object fit to a dataset
#' @param y Not used.
#' @param ... Not used.
#' @example man/example/bcm-plot_examples.R
#' @returns List: a list of two ggplot2 plots that are constructed using the information contained in an object of class bcm (Bayesian Clustering Model). A bcm object contains the clustering results from a clustering model that uses the TIP prior.
#' \item{trace_plot_posterior_number_of_clusters}{ggplot2 Plot: a plot of the posterior number of clusters (sampling iterations only) versus the corresponding sampling iteration number from the Gibbs sampler.}
#' \item{histogram_posterior_number_of_clusters}{ggplot2 Plot: a bar plot of the posterior number of clusters (sampling iterations only) from the Gibbs sampler.}
#' @exportMethod plot
setMethod("plot", signature(x="bcm",y="missing"), function(x,y,...){
  return(list(trace_plot_posterior_number_of_clusters = ggplot_number_of_clusters_trace(.posterior_number_of_clusters = x@posterior_number_of_clusters),
              histogram_posterior_number_of_clusters = ggplot_number_of_clusters_hist(.posterior_number_of_clusters = x@posterior_number_of_clusters)))
}
)
