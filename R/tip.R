# An R implementation of Bayesian clustering with the Table Invitation Prior (TIP)

# Global variables
globalVariables(c("%dopar%"))

# Read in some utility and plotting functions
# source("util.R")

#' @title Estimate The Number of Similar Subjects
#' @description Estimate the number of similar subjects \
#' using univariate multiple change point detection (i.e. binary segmentation).
#' @param .distance_matrix A symmetric n x n matrix of distances
#' @export
get_cpt_neighbors <- function(.distance_matrix){
  # --- A function used to obtain the nearest
  # neighbors for each subjects based on their
  # mutual distances ---

  # A vector to hold the number of subjects
  # most similar to subject i
  .num_neighbors <- vector()

  # Compute the number of subjects
  .n <- dim(.distance_matrix)[1]
  # Apply change-point detection to the set of
  # sorted distances corresponding to each subject i
  for(i in 1:.n){
    .num_neighbors[i] <- changepoint::cpt.mean(data = sort((.distance_matrix[i,])[-i],
                                                           decreasing = FALSE),
                                              method = "BinSeg", Q = floor(.n/2))@cpts[1]
  }
  return(.num_neighbors)
  # get_cpt_neighbors(.distance_matrix = matrix(abs(rnorm(100)),nrow=10,ncol=10))
}

#' @title Likelihood Models
#' @param .cluster_vector A vector of cluster assignments after the invitation step.
#' @param .i The subject index that the likelihood is computed for.
#' @param  .prior_estimates_for_likelihood A list of hyperparameters that are computed using
#' the data and are used in the likelihood function.
#' @param .likelihood_model The name of the likelihood model being used. Example: "NIW".
log_likelihood_fn <- function(.cluster_vector, .i, .prior_estimates_for_likelihood, .likelihood_model){
  if(toupper(.likelihood_model) == "NONE"){
    return(0)
  }else if(toupper(.likelihood_model == "NIW")){
    # Extract the prior parameters
    .data <- data.matrix(.prior_estimates_for_likelihood$.data)
    .lambda_0 <- as.numeric(.prior_estimates_for_likelihood$.lambda_0)
    .nu_0 <- as.numeric(.prior_estimates_for_likelihood$.nu_0)
    .mu_0 <- as.numeric(.prior_estimates_for_likelihood$.mu_0)
    .Psi_0 <- data.matrix(.prior_estimates_for_likelihood$.Psi_0)

    # number of variables for each subject
    .num_variables <- dim(.data)[2]

    # .n = number of subjects = sample size = number of rows in the dataset
    .n <- dim(.data)[1]

    # Extract the vector under consideration
    .yi <- .data[.i,]

    # Save the log-likelihood for each cluster
    .log_likelihood_vector <- vector()

    # Compute the number of clusters K
    if(length(.cluster_vector) == 0){
      return(0)
    }else{
      max_K <- max(.cluster_vector)
    }

    for(k in 1:max_K){
      # Find the subjects in the same cluster as subject .i (i.e. yi)
      .cluster_indices <- which(.cluster_vector == k)

      # Joint Prior: (mu, Sigma) ~ NIW(.mu_0, .lambda_0, .Psi_0, .nu_0)
      # Sampling Density: y_i|mu,Sigma ~ Np(mu,Sigma)
      # Joint Posterior: NIW(.mu_n, .lambda_n, Psi_n, .nu_n)

      # If the number of subjects in the same cluster as yi (including yi) is >= 3
      if(length(.cluster_indices) >= 3){
        .cluster_indices <- .cluster_indices[which(.cluster_indices != .i)]
      }else{
        # If the number of subjects in the same cluster as yi (including yi) is in {0,1,2},
        # then just set the cluster indices to 1, 2, ... .i - 1, .i + 1, .i + 2, ... , .n
        # since we need to have at least 2 subjects in each cluster (excluding yi) in order to
        # compute a log-likelihood
        .cluster_indices <- 1:.n
        .cluster_indices <- .cluster_indices[which(.cluster_indices != .i)]
      }

      # Isolate the subjects in the cluster
      .cluster_data <- .data[.cluster_indices,]

      # Compute the cluster size
      .n_k <- dim(.cluster_data)[1]

      # .ybar = column means where .ybar in R^p
      .ybar <- colMeans(.cluster_data)

      # .mu_n = (.lambda_0*.mu_0 + .n_k*.ybar)/(.lambda_0 + .n_k) where .mu_n in R^p
      .mu_n <- (.lambda_0*.mu_0 + .n_k*.ybar)/(.lambda_0 + .n_k)

      # .lambda_n = .lambda_0 + .n_k where .lambda_n in R^+
      .lambda_n = .lambda_0 + .n_k

      # .nu_n = .nu_0 + .n_k where .nu_n > p + 1
      .nu_n = .nu_0 + .n_k

      # .Psi_0 in R^{p x p} is the INVERSE scale matrix and is positive definite
      # Psi_n = .Psi_0 + S + (.lambda_0*.n_k)/(.lambda_0 + .n_k)(.ybar - .mu_0)(.ybar-.mu_0)^T
      # where S = sum_{i=1}^.n_k (y_i - .ybar)(y_i - .ybar)^T = Sigma*(.n_k-1) = cov(Y)*(.n_k-1), so we have
      # Psi_n = .Psi_0 + cov(.data)*(.n_k-1) + (.lambda_0*.n_k)/(.lambda_0 + .n_k)(.ybar - .mu_0)(.ybar-.mu_0)^T
      Psi_n = .Psi_0 + cov(.cluster_data)*(.n_k-1) + ((.lambda_0*.n_k)/(.lambda_0 + .n_k))*(.ybar - .mu_0)%*%t((.ybar-.mu_0))

      # Compute the scale matrix using the updated inverse scale matrix
      # If .S_temp is not symmetric, then round to make it symmetric
      # Note: symmetry of .S_temp is required for the Normal Inverse Wishart functions
      # Note: ---> thus we round to 5 digits for each matrix element
      .S_temp <- Psi_n #round(solve(Psi_n),digits = 5)

      # Draw mu and Sigma from their joint posterior distribution
      .draw <- LaplacesDemon::rnorminvwishart(n = 1,
                                              mu0 = as.numeric(.mu_n),
                                              lambda = as.numeric(.lambda_n),
                                              S = .S_temp,
                                              nu = as.numeric(.nu_n))

      # Add a small number to the diagonal until .S_temp is invertible
      # If .S_temp is already invertible, then do nothing
      .S_temp <- make_invertible(.matrix = .S_temp, .tolerance = 0.001)

      # Compute the Normal-Inverse-Wishart log-likelihood
      .log_likelihood_vector[k] <- LaplacesDemon::dnorminvwishart(mu = .yi,
                                                                  mu0 = .draw$mu,
                                                                  Sigma = .draw$Sigma,
                                                                  S = .S_temp,
                                                                  lambda = as.numeric(.lambda_n),
                                                                  nu = as.numeric(.nu_n),
                                                                  log = TRUE)

  }
    return(.log_likelihood_vector)
  }else{
      stop("Choose a valid likelihood function. Options are \"NIW\" and \"NONE\".")
    }
}

#' @title Table Invitation Prior Probability
#' @description Compute the prior probability that a subject belongs to a cluster.
#' @param .i The subject index (i.e. row index in a matrix for vector-variate data).
#' @param .similarity_matrix The matrix of similarity values.
#' @param .current_assignments The posterior assignments after the invitation step.
#' @param .num_clusters The number of clusters after the invitation step.
prob_tip_i <- function(.i, .similarity_matrix, .current_assignments, .num_clusters){
  # --- A function to compute the conditional posterior probability of a
  # subject joining a cluster (table) ---

  # Define a vector to hold the probability for subject i to be in cluster k = 1, 2, ..., <.num_clusters>
  .prob_k_vector <- vector()

  # Compute sum lambda(i,j) for j s.t. cluster[j] == k
  for(.k in 1:.num_clusters){
    .prob_k_vector[.k] <- sum(.similarity_matrix[.i, which(.current_assignments == .k)])
  }

  # Return the vector of probabilities
  return(.prob_k_vector)
}

#' @title Compute the Set of Similar Subjects
#' @description Find the <.num_candidates> subjects that are most similar to subject .i.
#' @param .i The subject index (i.e. row index in a matrix for vector-variate data).
#' @param .similarity_matrix The matrix of similarity values.
#' @param .num_candidates The number of similar subjects extracted.
get_candidates <- function(.i, .similarity_matrix, .num_candidates){
  # --- A function to return the <num_candidates> indices corresponding to
  # the subjects that are most similar to subject .i ---
  # Note: start at 2 since 1 is always the candidate itself
  return(order(.similarity_matrix[.i,], decreasing = TRUE)[2:(.num_candidates + 1)])
}

#' @title Bayesian Clustering with the Table Invitation Prior
#' @description Bayesian clustering with the Table Invitation Prior (TIP) and optional likelihood functions.
#' @param .data A list of vectors, matrices, or tensors that the analyst wishes to cluster.
#' @param .burn The number of burn-in iterations in the Gibbs sampler.
#' @param .samples The number of sampling iterations in the Gibbs sampler.
#' @param .similarity_matrix An n x x marix of simlarity values.
#' @param .init_num_neighbors A vector of integers corresponding to the estimate of the number of \
#'  number of subjects that are similar to the (i)th subject.
#'  @param .likelihood_model The name of the likelihood model used to compute the posterior probabilities.
#'  @param .subject_names An optional vector of names for the individual subjects.
#'  @param .num_cores The number of cores to use.
#'  @param .tolerance A parameter used to ensure matrices are invertible. A small number is iteratively \
#'  added to a matrix diagonal (if necessary) until the matrix is invertible.
#' @export
tip <- function(.data,
                .burn,
                .samples,
                .similarity_matrix,
                .init_num_neighbors,
                .likelihood_model = "NONE",
                .subject_names = vector(),
                .num_cores = 1,
                .tolerance = 0.001){

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
  }

  # Compute the total number of subjects
  .n <- dim(.similarity_matrix)[1]

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
  print(paste("Bayesian Clustering: Table Invitation Prior Gibbs Sampler"))
  print(paste("burn-in: ", .burn, sep = ""))
  print(paste("samples: ", .samples, sep = ""))
  print(paste("Likelihood Model: ", toupper(.likelihood_model), sep = ""))

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

    # Find the <p*.n> most similar neighbors to the random candidate
    .candidate_indices <- get_candidates(.i = .rand_init_candidate,
                                         .similarity_matrix = .similarity_matrix,
                                         .num_candidates = rpois(n = 1, lambda = .num_candidates))

    # The random candidate and its <p*.n> most similar subjects sit at the new table
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

      # Compute the conditional probabilities in parallel
      .posterior_assignment_temp <- foreach::foreach(.i = 1:.n) %dopar% {
        # Compute the log-prior for subject i for each cluster in the modified cluster vector
        .posterior_vector_i <- log(prob_tip_i(.i = .i,
                                              .similarity_matrix = .similarity_matrix,
                                              .num_clusters = .num_clusters,
                                              .current_assignments = .temp_cluster) + 1e-100)

        # Add the log-likelihood for subject i to the log prior
        .posterior_vector_i = .posterior_vector_i + log_likelihood_fn(.cluster_vector = .temp_cluster,
                                                                      .i = .i,
                                                                      .prior_estimates_for_likelihood = .prior_estimates_for_likelihood)

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
                                                                      .prior_estimates_for_likelihood = .prior_estimates_for_likelihood)
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

  # Return the necessary information to the analyst
  return(list(n = .n,
              burn = .burn,
              samples = .samples,
              posterior_assignments = .posterior_assignments[(.burn + 1):(.burn + .samples),],
              posterior_similarity_matrix = .posterior_similarity_matrix/.samples,
              posterior_number_of_clusters = .num_cluster_vector[(.burn + 1):(.burn + .samples)],
              prior_name = "TIP"))
}
