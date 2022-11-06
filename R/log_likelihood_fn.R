#' @title Likelihood Models
#' @param .cluster_vector Vector of positive integers. A vector of cluster assignments after the invitation step.
#' @param .i Positive integer. The integer subject index that the likelihood is computed for.
#' @param  .prior_estimates_for_likelihood List. A list of hyperparameters that are computed using
#' the data and are used in the likelihood function.
#' @param .likelihood_model Character. The character corresponding to the likelihood model being used. Options: "CONSTANT" (i.e. likelihood has no role in clustering), "NIW" for Normal-Inverse-Wishart likelihood model, "MNIW" for Matrix Normal-Inverse-Wishart likelihood model.
#' @importFrom stats cov rpois
#' @noRd
log_likelihood_fn <- function(.cluster_vector, .i, .prior_estimates_for_likelihood, .likelihood_model){
  if(toupper(.likelihood_model) == "CONSTANT"){
    return(0)
  }else if(toupper(.likelihood_model) == "NIW"){
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

      # .ybar = column means where .ybar in R^.p
      .ybar <- colMeans(.cluster_data)

      # .mu_n = (.lambda_0*.mu_0 + .n_k*.ybar)/(.lambda_0 + .n_k) where .mu_n in R^p
      .mu_n <- (.lambda_0*.mu_0 + .n_k*.ybar)/(.lambda_0 + .n_k)

      # .lambda_n = .lambda_0 + .n_k where .lambda_n in R^+
      .lambda_n = .lambda_0 + .n_k

      # .nu_n = .nu_0 + .n_k where .nu_n > .p + 1
      .nu_n = .nu_0 + .n_k

      # .Psi_0 in R^{.p x .p} is the INVERSE scale matrix and is positive definite
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
      .S_temp <- make_invertible(.matrix = .S_temp, .step_size = 0.001)

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
  }else if(toupper(.likelihood_model) == "MNIW"){
    # Extract the prior parameters
    .data <- .prior_estimates_for_likelihood$.data
    .row_cov <- .prior_estimates_for_likelihood$.row_cov
    .nu_c0 <- .prior_estimates_for_likelihood$.nu_c0
    .nu_r0 <- .prior_estimates_for_likelihood$.nu_r0
    .n <- .prior_estimates_for_likelihood$.n
    .m <- .prior_estimates_for_likelihood$.m
    .p <- .prior_estimates_for_likelihood$.p

    # Extract the subject matrix under consideration
    .yi <- .data[[.i]]

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

      # Isolate the subjects in the cluster  excluding the (.i)th subject
      .cluster_data <- .data[.cluster_indices]

      # ******************************************************************************************
      # *** NOTE: subject .i is excluded, so the computations below do not depend on subject .i ***
      # *** That is, although notation such as Ybar_k is used, in reality, it is actually Ybar_{k,-i}. ***
      # ******************************************************************************************

      # Compute the cluster size
      .n_k <- length(.cluster_data)

      # .Ybar_k: the m x .p mean matrix for the kth cluster
      .Ybar_k <- Reduce("+", .cluster_data)/length(.cluster_indices)
      # .ybar <- colMeans(.cluster_data)

      # .Psi_ck: a matrix parameter of the prior of Sigma_ck
      .Psi_ck <- cov(.Ybar_k)

      if(.p > 1){
        # Sigma_rk: covariance of the rows for the observations in the kth cluster
        .Sigma_rk <- cov(t(.Ybar_k))
      }else{
        # Special case: if .p = 1, then just use the overall row covariance matrix
        # In this case we can't compute cov(Ybar_k) since there is only 1 row
        .Sigma_rk <- .row_cov
      }

      # If .Sigma_rk is not invertible, then use the ridge trick to make it invertible;
      # that is, add the minimum amount of bias to make it invertible; otherwise do nothing.
      .Sigma_rk <- make_invertible_det(.Sigma_rk)

      # Compute the inverse of .Sigma_rk
      .Sigma_rk_inv = solve(.Sigma_rk)

      # .Lambda_k: mean of the Beta (i.e. a parameter referred to as "Beta") prior
      .Lambda_k = .Ybar_k

      # .Omega_k: row precision matrix of the Beta prior
      .Omega_k <- .Sigma_rk_inv#*.n_k

      # .Omega_k_hat: posterior row precision matrix of the
      # joint posterior distribution of .Beta_k and .Sigma_ck
      # (.Beta_k, .Sigma_ck) ~ MNIW(.Lambda_k_hat, .Omega_k_hat_inv, .Psi_ck_hat, .nu_ck)
      .Omega_k_hat <- .Sigma_rk_inv + .Omega_k

      # If .Omega_k_hat is not invertible, then use the ridge trick to make it invertible;
      # that is, add the minimum amount of bias to make it invertible.
      # Otherwise do nothing.

      .Omega_k_hat <- make_invertible_det(.Omega_k_hat)

      # .Omega_k_hat_inv: the inverse of .Omega_k_hat
      .Omega_k_hat_inv = solve(.Omega_k_hat)

      # .Lambda_k_hat: posterior mean of the joint posterior distribution of Beta_k and Sigma_ck
      # (.Beta_k, .Sigma_ck) ~ MNIW(.Lambda_k_hat, .Omega_k_hat_inv, .Psi_ck_hat, .nu_ck)
      .Lambda_k_hat = .Omega_k_hat_inv%*%(.Sigma_rk_inv%*%.Ybar_k + .Omega_k%*%.Lambda_k)

      # .Psi_ck_hat: posterior columns covariance of the joint posterior distribution of Beta_k and Sigma_ck
      # (.Beta_k, .Sigma_ck) ~ MNIW(.Lambda_k_hat, .Omega_k_hat_inv, .Psi_ck_hat, .nu_ck)
      .Psi_ck_hat = .Psi_ck + t(.Ybar_k) %*% .Sigma_rk_inv %*% .Ybar_k + t(.Lambda_k) %*% .Omega_k %*% .Lambda_k -
        t(.Lambda_k_hat) %*% .Omega_k_hat %*% .Lambda_k_hat

      # If Psi_ck_hat is not invertible, then use the ridge trick to make it invertible;
      # that is, add the minimum amount of bias to make it invertible.
      .Psi_ck_hat <- make_invertible_det(.Psi_ck_hat)

      # Compute .nu_ck
      .nu_ck = .nu_c0 + .m

      # Scale for identifiability
      .scale_factor_trace = .m/sum(diag(.Omega_k_hat))

      # Draw the posterior of Beta_k and Sigma_ck from the MNIW distribution
      .draw <- mniw::rmniw(n = 1,
                           Lambda = .Lambda_k_hat,
                           Omega = .scale_factor_trace*.Omega_k_hat,
                           Psi = .Psi_ck_hat,
                           nu = .nu_ck)

      # Extract the posterior draw of Beta_k
      .Beta_k_posterior = .draw$X

      # Extract the posterior draw of Sigma_ck
      .Sigma_ck_posterior = .draw$V

      # --- Compute the posterior of Sigma_rk ---
      # .Sigma_rk | Z_k ~ IW(.A_k + .Psi_rk, .nu_rk + .n)
      # Where .A_k = Z_k*Z_k^T and Z_k = (.Ybar_k - .mu_k)* .Sigma_ck^{-1/2}
      .svd_Sigma_ck_post <- svd(.Sigma_ck_posterior)

      # Compute Sigma_ck^{-1/2} (posterior Sigma_ck)
      if(.p==1){
        .Sigma_c_k_post_neg_half <- .svd_Sigma_ck_post$u%*%(.svd_Sigma_ck_post$d)^{-1/2}%*%.svd_Sigma_ck_post$v
      }else{
        .Sigma_c_k_post_neg_half <- .svd_Sigma_ck_post$u%*%diag((.svd_Sigma_ck_post$d)^{-1/2})%*%.svd_Sigma_ck_post$v
      }

      # Compute Z_k
      .Z_k = (.Ybar_k - .Beta_k_posterior)%*%.Sigma_c_k_post_neg_half

      # Compute A_k
      .A_k = .Z_k %*% t(.Z_k)

      # Initialize Psi_rk
      .Psi_rk = diag(.m)

      # Compute nu_rk
      .nu_rk = .nu_r0 + .p

      # Draw the posterior Sigma_rk | Z_k
      .Sigma_rk_posterior = LaplacesDemon::rinvwishart(nu = .nu_rk, S = .A_k + .Psi_rk)

      # Compute the likelihood on the log-scale assuming a Matrix Normal Distribution
      .log_likelihood_vector[k] <- LaplacesDemon::dmatrixnorm(X = .data[[.i]],
                                                              M = .Beta_k_posterior,
                                                              U = .n_k*.Sigma_rk_posterior,
                                                              V = .n_k*.Sigma_ck_posterior,
                                                              log = TRUE)

    }
    return(.log_likelihood_vector)
  }else{
    stop("Choose a valid likelihood function. Options are \"NIW\" and \"CONSTANT\".")
  }
}
