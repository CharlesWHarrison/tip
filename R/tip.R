#' @title Bayesian Clustering with the Table Invitation Prior
#' @description Bayesian clustering with the Table Invitation Prior (TIP) and optional likelihood functions.
#' @param .data Data frame (vectors comprise a row in a data frame; NIW only) or a list of matrices (MNIW only) that the analyst wishes to cluster. Note: if .likelihood_model = "CONSTANT", then the .data argument has no effect.
#' @param .burn Non-negative integer: the number of burn-in iterations in the Gibbs sampler.
#' @param .samples Positive integer: the number of sampling iterations in the Gibbs sampler.
#' @param .similarity_matrix Matrix: an n x n matrix of similarity values.
#' @param .init_num_neighbors Vector of positive integers: each (i)th positive integer corresponds to the estimate of the number of subjects that are similar to the (i)th subject.
#' @param .likelihood_model Character: the name of the likelihood model used to compute the posterior probabilities. Options: "NIW" (vectors; .data is a dataframe), "MNIW" (matrices; .data is a list of matrices), or "CONSTANT" (vector, matrices, and tensors; .data is an empty list)
#' @param .subject_names Vector of characters: an optional vector of names for the individual subjects. This is useful for the plotting function.
#' @param .num_cores Positive integer: the number of cores to use.
#' @param .step_size Positive numeric: A parameter used to ensure matrices are invertible. A small number is iteratively added to a matrix diagonal (if necessary) until the matrix is invertible.
#' @returns Object of class bcm: bcm denotes a "Bayesian Clustering Model" object that contains the results from a clustering model that uses the TIP prior.
#' \item{n}{Positive integer: the sample size or number of subjects.}
#' \item{burn}{Non-negative integer: the number of burn-in iterations in the Gibbs sampler.}
#' \item{samples}{Positive integer: the number of sampling iterations in the Gibbs sampler.}
#' \item{posterior_assignments}{List: a list of <\code{samples}> vectors where the ith element of each vector is the posterior cluster assignment for the ith subject.}
#' \item{posterior_similarity_matrix}{Matrix: an \code{n} x \code{n} matrix where the (i,j)th element is the posterior probability that subject i and subject j are in the same cluster.}
#' \item{posterior_number_of_clusters}{Vector of positive integers: a vector where the jth element is the number of clusters after posterior sampling (i.e., the posterior number of clusters).}
#' \item{prior_name}{Character: the name of the prior used.}
#' \item{likelihood_name}{Character: the name of the likelihood used.}
#' @example man/example/tip_examples.R
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom foreach %dopar%
#' @importFrom methods new
#' @export
tip <- function(.data = list(),
                .burn = 1000,
                .samples = 1000,
                .similarity_matrix,
                .init_num_neighbors,
                .likelihood_model = "CONSTANT",
                .subject_names = vector(),
                .num_cores = 1,
                .step_size = 0.001){

  # Some initial checks and informative error message.
  if(.burn < 0){
    stop(".burn must be an integer >= 0.")
  }else if(.samples < 0){
    stop(".samples must be an integer > 0.")
  }else if(.likelihood_model %in% c("NIW", "MNIW", "CONSTANT") == FALSE){
    stop(".likelihood_model must be one of the following characters
         'NIW', 'MNIW', or 'CONSTANT'.")
  }else if(.num_cores <= 0){
    stop(".num_cores must be an integer > 0 like 1, 2, 3, ....")
  }else if(.step_size <= 0){
    stop(".step_size must be a float > 0.")
  }else{1}

  # Compute the total number of subjects
  .n <- dim(.similarity_matrix)[1]

  if(toupper(.likelihood_model) == "NIW"){
    # NIW: Normal-Inverse-Wishart
    .Psi_0 <- cov(.data)*(dim(.data)[2]-1)
    .Psi_0 <- solve(.Psi_0)
    .Psi_0 <- (.Psi_0 + t(.Psi_0))/2
    .lambda_0 <- 1
    .nu_0 <- dim(.data)[1] + 1
    .mu_0 <- as.numeric(colMeans(.data))
    .prior_estimates_for_likelihood = list(.data = .data,
                                           .Psi_0 = .Psi_0,
                                           .lambda_0 = .lambda_0,
                                           .nu_0 = .nu_0,
                                           .mu_0 = .mu_0)
  }else if(toupper(.likelihood_model) == "MNIW"){

    # Average the Y values over all list elements
    # .Ybar is .m x .p and each Yi is .m x .p
    # This checks if all matrices are the same size
    .Ybar <- Reduce('+', .data)/length(.data)

    # Compute the number of rows and columns in the response matrix Y
    .dim_matrices = dim(.Ybar)

    # Let n = the number of response matrices Yi for i = 1, 2, ..., n
    .n = length(.data)

    # Let m = # of rows in each response matrix Yi = # rows in the X matrix
    .m = .dim_matrices[1]

    # Let .p = # of columns in each response matrix Yi
    .p = .dim_matrices[2]

    # When there is only a single variable, then .p = 1 and we have a set of
    # n vectors with dimension m x 1
    .row_cov = NA
    if(.p ==1 & .m == .n){
      .row_cov = cov(data.frame(do.call("rbind", .data)))
    }

    if(.p == 1 & .m != .n){
      .row_cov = diag(.m)
    }

    # .nu_c0 and .nu_r0 are the prior degrees of freedom
    # .nu_c0 and .nu_r0 are the prior degrees of freedom
    if(.n < .m & .n < .p){
      .nu_c0 <- .n + max(.m,.p)
      .nu_r0 <- .n + max(.m,.p)
    }else if(.n < .m & .n >= .p){
      .nu_c0 <- .n + max(.m,.p)
      .nu_r0 <- .n + max(.m,.p)
    }else if(.n >= .m & .n >= .p){
      .nu_c0 <- .n
      .nu_r0 <- .n
    }else if(.n >= .m & .n < .p){
      .nu_c0 <- .n + max(.m,.p)
      .nu_r0 <- .n + max(.m,.p)
    }else{
      stop("Failure to set degrees of freedom nu_c0 and nu_r0")
    }

    # Create a list of prior parameters
    .prior_estimates_for_likelihood = list(.data = .data,
                                           .row_cov = .row_cov,
                                           .p = .p,
                                           .m = .m,
                                           .nu_c0 = .nu_c0,
                                           .nu_r0 = .nu_r0,
                                           .n = .n)
  }

  # Store the cluster assignments
  .posterior_assignments <- list()

  # In the first iteration there is exactly 1 cluster
  .posterior_assignments[[1]] <- rep(1,.n)

  # Vector to store number of clusters
  .num_cluster_vector <- vector()
  .num_cluster_vector[1] <- length(table(.posterior_assignments))

  # Intialize the posterior similarity matrix
  .posterior_similarity_matrix <- matrix(0, nrow = .n, ncol = .n)

  # Create a progress bar
  .tip_cpt_pb = txtProgressBar(min = 2, max = .burn + .samples, style = 3)

  # Print message to the analyst
  message(paste("Bayesian Clustering: Table Invitation Prior Gibbs Sampler"))
  message(paste("burn-in: ", .burn, sep = ""))
  message(paste("samples: ", .samples, sep = ""))
  message(paste("Likelihood Model: ", toupper(.likelihood_model), sep = ""))

  # Iteration .t = 1, 2, ..., .burn + .samples gives <.samples> .samples from the posterior
  for(.t in 2:(.burn + .samples)){
    # Update the progress bar
    setTxtProgressBar(.tip_cpt_pb, .t)

    # Extract the cluster assignments from the previous iteration
    .temp_cluster <- .posterior_assignments[[.t-1]]

    # Compute the current number of clusters
    .num_clusters <- length(unique(.posterior_assignments[[.t-1]]))

    # Save the number of clusters
    .num_cluster_vector[.t] <- .num_clusters

    # Pick a random candidate for the new cluster (table)
    .rand_init_candidate <- sample(x = 1:.n, size = 1, replace = FALSE, prob = rep(1/.n,.n))

    # The number of candidates is based on the number of subjects to the left
    # of the first change-point
    .num_candidates <- .init_num_neighbors[.rand_init_candidate]

    # The number of candidates must be between two and .n - 1 in order to compute
    # the necessary inverse matrices, etc.
    if(.num_candidates < 2){.num_candidates <- 2}
    if(.num_candidates == .n){.num_candidates = .n - 1}

    # Find the <.p*.n> most similar neighbors to the random candidate
    .candidate_indices <- get_candidates(.i = .rand_init_candidate,
                                         .similarity_matrix = .similarity_matrix,
                                         .num_candidates = rpois(n = 1, lambda = .num_candidates))

    # The random candidate and its <.p*.n> most similar subjects sit at the new table
    .temp_cluster[c(.rand_init_candidate, .candidate_indices)] <- .num_clusters + 1

    # Ensure that the cluster assignments are contiguous: 1, 2, ..., K
    .temp_cluster <- recode(.temp_cluster)

    # Recompute the number of clusters now that a new cluster has been added
    .num_clusters <- length(table(.temp_cluster))

    # For each subject .t = 1, 2, ..., .n compute the posterior probability for each cluster
    # independently of the other subjects and sample the posterior assignment
    .posterior_assignment_temp <- vector()

    if(.num_cores > 1){
      # library(foreach)
      # Allocate the logical processors
      cl <- parallel::makeCluster(.num_cores)
      doParallel::registerDoParallel(cl)

      # Export the following functions to each core
      parallel::clusterExport(cl,list('prob_tip_i',
                                      'log_likelihood_fn',
                                      'make_invertible'),
                              envir=environment())

      # Export the following libraries to each core
      parallel::clusterEvalQ(cl, c(library(LaplacesDemon)))

      # Initial value for i so that the warning from roxygen2 disappears
      i <- 1

      # Compute the conditional probabilities in parallel
      .posterior_assignment_temp <- foreach::foreach(i = 1:.n) %dopar%{
        # Compute the log-prior for subject .i for each cluster in the modified cluster vector
        .posterior_vector_i <- log(prob_tip_i(.i = i,
                                              .similarity_matrix = .similarity_matrix,
                                              .num_clusters = .num_clusters,
                                              .current_assignments = .temp_cluster) + 1e-100)

        # Add the log-likelihood for subject i to the log prior
        .posterior_vector_i = .posterior_vector_i + log_likelihood_fn(.cluster_vector = .temp_cluster,
                                                                      .i = i,
                                                                      .prior_estimates_for_likelihood = .prior_estimates_for_likelihood,
                                                                      .likelihood_model = .likelihood_model)

        # Convert to a posterior probability
        .posterior_vector_i <- sapply(.posterior_vector_i, function(qq) exp(qq - max(.posterior_vector_i)))
        .posterior_vector_i <- sapply(.posterior_vector_i, function(qq) qq/sum(.posterior_vector_i))

        # Sample subject i's posterior cluster assignment from the posterior
        # Note: since this is the last line in the foreach loop, the result
        # is saved as an element in a list
        sample(x = 1:(.num_clusters), size = 1, replace = FALSE, prob = .posterior_vector_i)
      }
      # Deallocate the parallel resources
      parallel::stopCluster(cl)

      # Convert from list to a vector
      .posterior_assignment_temp <- unlist(.posterior_assignment_temp)

    }else{
      for(.i in 1:.n){
        # Compute the log-prior for subject i for each cluster in the modified cluster vector
        .posterior_vector_i <- log(prob_tip_i(.i = .i,
                                              .similarity_matrix = .similarity_matrix,
                                              .num_clusters = .num_clusters,
                                              .current_assignments = .temp_cluster) + 1e-100)

        # Add the log-likelihood for subject i to the log prior
        .posterior_vector_i = .posterior_vector_i + log_likelihood_fn(.cluster_vector = .temp_cluster,
                                                                      .i = .i,
                                                                      .prior_estimates_for_likelihood = .prior_estimates_for_likelihood,
                                                                      .likelihood_model = .likelihood_model)
        # Convert to a posterior probability
        .posterior_vector_i <- sapply(.posterior_vector_i, function(qq) exp(qq - max(.posterior_vector_i)))
        .posterior_vector_i <- sapply(.posterior_vector_i, function(qq) qq/sum(.posterior_vector_i))

        # Sample subject i's posterior cluster assignment from the posterior
        .posterior_assignment_temp <- c(.posterior_assignment_temp,
                                        sample(x = 1:(.num_clusters),
                                               size = 1,
                                               replace = FALSE,
                                               prob = .posterior_vector_i))
      }

    }
    # Recode the posterior probability so that the cluster assignments
    # have the form 1, 2, ..., K (i.e. contiguous values)
    .posterior_assignments[[.t]] <- recode(.posterior_assignments = .posterior_assignment_temp)

    # Update the posterior similarity matrix
    # *** NOTE: do not divide by the number of iterations <.t>
    if(.t > .burn){
      .posterior_similarity_matrix = .posterior_similarity_matrix + get_proximity_matrix(.assignments = .posterior_assignments[[.t]])
    }
  }

  # Close the progress bar
  close(.tip_cpt_pb)

  # Put the posterior cluster assignments in a data frame
  .posterior_assignments <- do.call(rbind.data.frame, .posterior_assignments)
  colnames(.posterior_assignments) <- ifelse(length(.subject_names) == 0, paste("Subject", 1:.n), .subject_names)

  return(new("bcm",
             n = .n,
             burn = .burn,
             samples = .samples,
             posterior_assignments = .posterior_assignments[(.burn + 1):(.burn + .samples),],
             posterior_similarity_matrix = .posterior_similarity_matrix/.samples,
             posterior_number_of_clusters = .num_cluster_vector[(.burn + 1):(.burn + .samples)],
             prior_name = "TIP",
             likelihood_name = .likelihood_model))
}
