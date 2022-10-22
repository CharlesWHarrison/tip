#' @title Bayesian Clustering Model (bcm) S4 class.
#' @description An S4 class to store the results of the Gibbs sampler.
#' @slot n Positive integer; the sample size (i.e. number of subjects).
#' @slot burn Non-negative integer; the number of burn-in iterations in the Gibbs sampler.
#' @slot samples Positive integer; the number of sampling iterations in the Gibbs sampler.
#' @slot posterior_assignments List; a list of vectors of cluster assignments (positive integers) for each sampling iteration in the Gibbs sampler.
#' @slot posterior_similarity_matrix Matrix; a matrix where the (i,j)th element is the posterior probability that subject i and subject j are in the same cluster.
#' @slot posterior_number_of_clusters Vector of positive integers; each vector element is the number of clusters after posterior sampling for each sampling iteration in the Gibbs sampler.
#' @slot prior_name Character. The name of the prior used.
#' @exportClass bcm
setClass("bcm",
         slots=list(n = "numeric",
                    burn = "numeric",
                    samples = "numeric",
                    posterior_assignments = "data.frame",
                    posterior_similarity_matrix = "matrix",
                    posterior_number_of_clusters = "numeric",
                    prior_name = "character"))

#' @title A function to return plots from a Bayesian Clustering Model (bcm) object
#' @rdname plot.bcm
#' @aliases plot
#' @param x A Bayesian Clustering Model (bcm) object
#' @param y Not used.
#' @param ... Not used.
#' @exportMethod plot
setMethod("plot", signature(x="bcm",y="missing"), function(x,y,...){
  return(list(trace_plot_posterior_number_of_clusters = ggplot_number_of_clusters_trace(.posterior_number_of_clusters = x@posterior_number_of_clusters),
              histogram_posterior_number_of_clusters = ggplot_number_of_clusters_hist(.posterior_number_of_clusters = x@posterior_number_of_clusters)))
}
)
