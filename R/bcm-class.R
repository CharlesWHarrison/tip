#' @title Bayesian Clustering Model (bcm) S4 class.
#' @description An S4 class to store the results of the Gibbs sampler.
#' @slot n Positive integer: the sample size (i.e., the number of subjects).
#' @slot burn Non-negative integer: the number of burn-in iterations in the Gibbs sampler.
#' @slot samples Positive integer: the number of sampling iterations in the Gibbs sampler.
#' @slot posterior_assignments List of vectors of positive integers: a list of vectors of cluster assignments (i.e., positive integers) for each sampling iteration in the Gibbs sampler.
#' @slot posterior_similarity_matrix Matrix: a matrix where the (i,j)th element is the posterior probability that subject i and subject j belong to the same cluster.
#' @slot posterior_number_of_clusters Vector of positive integers: each vector element is the number of clusters after posterior sampling for each sampling iteration in the Gibbs sampler.
#' @slot prior_name Character: the name of the prior used.
#' @exportClass bcm
setClass("bcm",
         slots=list(n = "numeric",
                    burn = "numeric",
                    samples = "numeric",
                    posterior_assignments = "data.frame",
                    posterior_similarity_matrix = "matrix",
                    posterior_number_of_clusters = "numeric",
                    prior_name = "character",
                    likelihood_name = "character"))


