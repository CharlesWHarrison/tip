# Clustering vectors, matrices, and tensors using the Table Invitation Prior (TIP) in R   
## Installation 
``` 
# Install from CRAN
install.packages("tip") 
```

## Paper Citation
Charles W. Harrison, Qing He, and Hsin-Hsiung Huang. “Clustering Gene Expressions Using the Table Invitation Prior”. In: Genes 13.11 (2022). issn: 2073-4425. doi: 10.3390/genes13112036. url: https://www.mdpi.com/2073-4425/13/11/2036.

```text
@Article{genes13112036,
AUTHOR = {Harrison, Charles W. and He, Qing and Huang, Hsin-Hsiung},
TITLE = {Clustering Gene Expressions Using the Table Invitation Prior},
JOURNAL = {Genes},
VOLUME = {13},
YEAR = {2022},
NUMBER = {11},
ARTICLE-NUMBER = {2036},
URL = {https://www.mdpi.com/2073-4425/13/11/2036},
ISSN = {2073-4425},
DOI = {10.3390/genes13112036}
}
```

## R Package Citation
```text
@Manual{tip2022,
    title = {tip: Bayesian Clustering Using the Table Invitation Prior (TIP)},
    author = {Charles W. Harrison and Qing He and Hsin-Hsiung Huang},
    year = {2022},
    note = {R package version 0.1.0},
    url = {https://CRAN.R-project.org/package=tip},
  }
```

## Introduction
This R library provides a Gibbs sampler for Bayesian clustering models that utilize the Table Invitation Prior (TIP) introduced by Harrison, He, and Huang (2022). TIP utilizes pairwise distance and pairwise similarity information between the observed data (i.e., subjects). The term ''subject'' is used to refer to an individual vector, matrix, or higher-order tensors. 
1. **Hypothetical Vector-variate subject example**: doctors measure the 5 vital signs of 19 people and thus there are 19 subjects that are each described by a ``5 x 1 `` vector. Each individual person is considered to be a subject, so there are 19 total subjects. 
2. **Hypothetical Matrix-variate subject example**: a single X-ray is taken for 57 adults, and each X-ray image is stored as a ``512 x 512`` matrix where each value in each matrix varies between zero and one (i.e., a grayscale image). Each X-ray is considered as an individual subject, so there are 57 subjects.
3. **Hypothetical Tensor-variate subject example**: 23 adults have an fMRI taken. Each fMRI image corresponds to a 3-way tensor, each of the 23 3-way tensors correspond to an individual subject, so there are 23 subjects.

In Bayesian clustering, the goal is sample from the posterior distribution which is given by 
$$P(\mathbf{c}|\mathbf{x}) \propto P(\mathbf{X}|\mathbf{c})P(\mathbf{c})$$ 
where $P(\mathbf{X}|\mathbf{c})$ is the likelihood function and $P(\mathbf{c})$ is the prior distribution. In this case, $P(\mathbf{c})$ refers to the Table Invitation Prior (TIP). TIP is flexible and may be incorporated with a variety of different likelihood functions designed for different types of data. The current options for the likliehood function are the following:

1. The ``.likelihood_model = "CONSTANT"`` is the fastest option and can be used for vectors, matrices, and higher-order tensors (i.e., ``.data`` is not used). The "CONSTANT" option returns a constant likelihood function value regardless of the observed data so that likelihood function has no role in the clustering. The "CONSTANT" likelihood option may be used for vector-variate datasets (e.g. the Iris dataset, US Arrests dataset, etc.), matrix-variate datasets (e.g. data pertaining to electroencephalograms (EEGs), grayscale images, etc.), and higher-order tensor-variate datasets (i.e., videos, colored-pictures, etc.). 

2. The ``.likelihood_model = "NIW"`` option can be used for vector-datasets only (i.e., ``.data`` is a ``data.frame``). The "NIW" option uses a "Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of vector-variate data include the Iris dataset, US Arrests dataset, and so on.  

3. The ``.likelihood_model = "MNIW"`` option can be used for matrix-variate data (i.e., ``.data`` is a ``list`` of matrices). The "MNIW" option uses a "Matrix-Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of matrix-variate data include EEG data, grayscale images (i.e., X-rays or black and white photographs), graph data, and so on. Note that "MNIW" may be used for vector-variate data by passing a list of matrices that have either 1 row and multiple columns or matrices that have one column and multiple rows.

## Clustering the Iris Dataset (i.e., vectors) with a Normal-Inverse-Wishart (NIW) likelihood and a TIP prior
```
  # Import the tip library
  library(tip)

  # Import the iris dataset
  data(iris)

  # The first 4 columns are the data whereas
  # the 5th column refers to the true labels
  X <- data.matrix(iris[,c("Sepal.Length",
                           "Sepal.Width",
                           "Petal.Length",
                           "Petal.Width")])

  # Extract the true labels (optional)
  # True labels are only necessary for constructing network
  # graphs that incorporate the true labels; this is often
  # for research.
  true_labels <- iris[,"Species"]

  # Compute the distance matrix
  distance_matrix <- data.matrix(dist(X))

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMMENDATION: burn >= 1000
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMMENDATION: samples >= 1000
  samples <- 1000

  # Set the subject names
  names_subjects <- paste(1:dim(iris)[1])

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = data.matrix(X),
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "NIW",
              .subject_names = names_subjects,
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)

  # View the posterior distribution of the number of clusters
  tip_plots$histogram_posterior_number_of_clusters

  # View the trace plot with respect to the posterior number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
  cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

  # If the true labels are available, then show the cluster result via a contingency table
  table(data.frame(true_label = true_labels,
                   cluster_assignment = cluster_assignments))

  # Create the one component graph with minimum entropy
  partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                               .num_components = 1,
                                               .step_size = 0.001)

  # Associate class labels and colors for the plot
  class_palette_colors <- c("setosa" = "blue",
                            "versicolor" = 'green',
                            "virginica" = "orange")

  # Associate class labels and shapes for the plot
  class_palette_shapes <- c("setosa" = 19,
                            "versicolor" = 18,
                            "virginica" = 17)

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = NA,
                      .subject_class_names = true_labels,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = FALSE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels. Also, suppress
  # the subject labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = FALSE)
```
## Clustering the US Arrests Dataset (i.e., vectors) with a Normal-Inverse-Wishart (NIW) likelihood and a TIP prior
```
 # Import the TIP library
  library(tip)

  # Import the US Arrests dataset
  data(USArrests)

  # Convert the data to a matrix
  X <- data.matrix(USArrests)

  # Compute the distance matrix
  distance_matrix <- data.matrix(dist(X))

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMENDATION: *** burn >= 1000 ***
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMENDATION: *** samples >= 1000 ***
  samples <- 1000

  # Extract the state names
  names_subjects <- rownames(USArrests)

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = data.matrix(X),
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "NIW",
              .subject_names = names_subjects,
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)

  # View the posterior distribution of the number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # View the trace plot with respect to the posterior number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
  cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

  # Create a list where each element stores the cluster assignments
  cluster_assignment_list <- list()
  for(k in 1:length(unique(cluster_assignments))){
    cluster_assignment_list[[k]] <- names_subjects[which(cluster_assignments == k)]
  }
  cluster_assignment_list

  # Create the one component graph with minimum entropy
  partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                               .num_components = 1,
                                               .step_size = 0.001)

  # View the state names
  # names_subjects

  # Create a vector of true region labels to see if there is a pattern
  true_region <- c("Southeast", "West", "Southwest", "Southeast", "West", "West",
                   "Northeast", "Northeast", "Southeast", "Southeast", "West", "West",
                   "Midwest", "Midwest", "Midwest", "Midwest", "Southeast", "Southeast",
                   "Northeast", "Northeast", "Northeast", "Midwest", "Midwest", "Southeast",
                   "Midwest", "West", "Midwest", "West", "Northeast", "Northeast",
                   "Southwest", "Northeast", "Southeast", "Midwest", "Midwest", "Southwest",
                   "West", "Northeast", "Northeast", "Southeast", "Midwest", "Southeast",
                   "Southwest", "West", "Northeast", "Southeast", "West", "Southeast",
                   "Midwest", "West")

  names_subjects

  # Associate class labels and colors for the plot
  class_palette_colors <- c("Northeast" = "blue",
                            "Southeast" = "red",
                            "Midwest" = "purple",
                            "Southwest" = "orange",
                            "West" = "green")

  # Associate class labels and shapes for the plot
  class_palette_shapes <- c("Northeast" = 15,
                            "Southeast" = 16,
                            "Midwest" = 17,
                            "Southwest" = 18,
                            "West" = 19)

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .subject_class_names = true_region,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = TRUE)


  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  # Remove the subject names with .add_node_labels = FALSE
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .subject_class_names = true_region,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = FALSE)

  # Construct a network plot without class labels
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # Construct a network plot without class labels. Also, suppress
  # the subject labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = FALSE)
```
## Clustering gene expression data (i.e., vectors) with a Normal-Inverse-Wishart (NIW) likelihood and a TIP prior
```
 # Import the TIP library
  library(tip)

  # ----- Dataset information -----
  # The data were accessed from the UCI machine learning repository
  # Original link: https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq
  # Source: Samuele Fiorini, samuele.fiorini '@' dibris.unige.it, University of Genoa, redistributed under Creative Commons license (http://creativecommons.org/licenses/by/3.0/legalcode) from https://www.synapse.org/#!Synapse:syn4301332.
  # Data Set Information: Samples (instances) are stored row-wise. Variables (attributes) of each sample are RNA-Seq gene expression levels measured by illumina HiSeq platform.
  # Relevant Papers: Weinstein, John N., et al. 'The cancer genome atlas pan-cancer analysis project.' Nature genetics 45.10 (2013): 1113-1120.
  # -------------------------------

  # Import the data (see the provided link above)
  X <- read.csv("data.csv")
  true_labels <- read.csv("labels.csv")

  # Extract the true indices
  subject_names <- true_labels$X

  # Extract the true classes
  true_labels <- true_labels$Class

  # Convert the dataset into a matrix
  X <- data.matrix(X)

  ##### BEGIN: Apply PCA to the dataset #####

  # Step 1: perform Prinicpal Components Analysis (PCA) on the dataset
  pca1 <- prcomp(X)

  # Step 2: compute summary information
  summary1 <- summary(pca1)

  # Step 3: plot the cumulative percentage of the variance
  # explained against the number of principal components
  tip::ggplot_line_point(.x = 1:length(summary1$importance['Cumulative Proportion',]),
                         .y = summary1$importance['Cumulative Proportion',],
                         .xlab = "Number of Principal Components",
                         .ylab = "Cumulative Percentage of the Variance Explained")

  # The number of principal components chosen is 7, and
  # the 7 principal components explain roughly 80% of
  # the variance
  num_principal_components <- which(summary1$importance['Cumulative Proportion',] <= 0.8)

  # The clustering is applied to the principal component dataset
  X <- pca1$x[,as.numeric(num_principal_components)]
  ##### END: Apply PCA to the dataset #####

  # Compute the distance matrix
  distance_matrix <- data.matrix(dist(X))

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMENDATION: *** burn >= 1000 ***
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMENDATION: *** samples >= 1000 ***
  samples <- 1000

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = data.matrix(X),
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "NIW",
              .subject_names = c(),
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)

  # View the posterior distribution of the number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # View the trace plot with respect to the posterior number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
  cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

  # Create a list where each element stores the cluster assignments
  cluster_assignment_list <- list()
  for(k in 1:length(unique(cluster_assignments))){
    cluster_assignment_list[[k]] <- true_labels[which(cluster_assignments == k)]
  }
  cluster_assignment_list

  # Create the one component graph with minimum entropy
  partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                               .num_components = 1,
                                               .step_size = 0.001)

  # Associate class labels and shapes for the plot
  class_palette_shapes <- c("PRAD" = 19,
                            "BRCA" = 18,
                            "KIRC" = 17,
                            "LUAD" = 16,
                            "COAD" = 15)

  # Associate class labels and colors for the plot
  class_palette_colors <- c("PRAD" = "blue",
                            "BRCA" = "red",
                            "KIRC" = "black",
                            "LUAD" = "green",
                            "COAD" = "orange")

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = subject_names,
                      .subject_class_names = true_labels,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  # Remove the subject names with .add_node_labels = FALSE
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = subject_names,
                      .subject_class_names = true_labels,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = FALSE)

  # Construct a network plot without class labels but with subject labels
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = subject_names,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # Construct a network plot without class labels and subject labels
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = subject_names,
                      .node_size = 2,
                      .add_node_labels = FALSE)
```

## Clustering Matrix-variate Data with a Matrix-Normal-Inverse-Wishart (MNIW) likelihood and a TIP prior
```
   # Import the tip library 
  library(tip)

  # A function to generate random matrices from a matrix normal distribution
  random_mat_normal <- function(mu, num_rows, num_cols){
    LaplacesDemon::rmatrixnorm(M = matrix(mu,
                                          nrow = num_rows,
                                          ncol = num_cols),
                               U = diag(num_rows),
                               V = diag(num_cols))
  }

  # Generate 3 clusters of matrices
  p <- 5
  m <- 3
  c1 <- lapply(1:20, function(x) random_mat_normal(mu = 0, num_rows = m, num_cols = p))
  c2 <- lapply(1:25, function(x) random_mat_normal(mu = 5, num_rows = m, num_cols = p))
  c3 <- lapply(1:30, function(x) random_mat_normal(mu = -5, num_rows = m, num_cols = p))

  # Put all the data into a list
  data_list <- c(c1,c2,c3)

  # Create a vector of true labels. True labels are only necessary
  # for constructing network graphs that incorporate the true labels;
  # this is often useful for research.
  true_labels <- c(rep("Cluster 1", 20),
                   rep("Cluster 2", 25),
                   rep("Cluster 3", 30))

  distance_matrix <- matrix(NA,
                            nrow = length(true_labels),
                            ncol = length(true_labels))
  # Distance matrix
  for(i in 1:length(true_labels)){
    for(j in i:length(true_labels)){
      distance_matrix[i,j] <- SMFilter::FDist2(mX = data_list[[i]],
                                               mY = data_list[[j]])
      distance_matrix[j,i] <- distance_matrix[i,j]
    }
  }

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMMENDATION: burn >= 1000
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMMENDATION: samples >= 1000
  samples <- 1000

  # Set the subject names
  names_subjects <- paste(1:dim(distance_matrix)[1])

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = data_list,
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "MNIW",
              .subject_names = names_subjects,
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)

  # View the posterior distribution of the number of clusters
  tip_plots$histogram_posterior_number_of_clusters

  # View the trace plot with respect to the posterior number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
  cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

  # If the true labels are available, then show the cluster result via a contigency table
  table(data.frame(true_label = true_labels,
                   cluster_assignment = cluster_assignments))

  # Create the one component graph with minimum entropy
  partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                               .num_components = 1,
                                               .step_size = 0.001)

  # Associate class labels and colors for the plot
  class_palette_colors <- c("Cluster 1" = "blue",
                            "Cluster 2" = 'green',
                            "Cluster 3" = "red")

  # Associate class labels and shapes for the plot
  class_palette_shapes <- c("Cluster 1" = 19,
                            "Cluster 2" = 18,
                            "Cluster 3" = 17)

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = NA,
                      .subject_class_names = true_labels,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = FALSE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels. Also, suppress the
  # subject labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = FALSE)
```

## Clustering Tensor-variate Data with CONSTANT likelihood and a TIP prior
```
# Import the tip library
  library(tip)

  # ################## NOTE ##################
  # Order 3 Tensor dimension: c(d1,d2,d3)
  # d1: number of rows
  # d2: number of columns
  # d3: number of slices
  ############################################

  # Set a random seed for reproducibility
  set.seed(007)

  # A function to generate an order-3 tensor
  generate_gaussian_tensor <- function(.tensor_dimension, .mean = 0, .sd = 1){
    array(data = c(rnorm(n = prod(.tensor_dimension),
                         mean = .mean,
                         sd = .sd)),
          dim = .tensor_dimension)
  }

  # Define the tensor dimension
  tensor_dimension <- c(256,256,3)

  # Generate clusters of tensors
  c1 <- lapply(1:30, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                          .mean = 0,
                                                          .sd = 1))


  # Generate clusters of tensors
  c2 <- lapply(1:40, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                          .mean = 5,
                                                          .sd = 1))


  # Generate clusters of tensors
  c3 <- lapply(1:50, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                          .mean = -5,
                                                          .sd = 1))

  # Make a list of tensors
  X <- c(c1, c2, c3)

  # Compute the number of subjects for each cluster
  n1 <- length(c1)
  n2 <- length(c2)
  n3 <- length(c3)

  # Create a vector of true labels. True labels are only necessary
  # for constructing network graphs that incorporate the true labels;
  # this is often useful for research.
  true_labels <- c(rep("Cluster 1", n1),
                   rep("Cluster 2", n2),
                   rep("Cluster 3", n3))

  # Compute the total number of subjects
  n <- length(X)

  # Construct the distance matrix
  distance_matrix <- matrix(data = NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in i:n){
      distance_matrix[i,j] <- sqrt(sum((X[[i]] - X[[j]])^2))
      distance_matrix[j,i] <- distance_matrix[i,j]
    }
  }

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMMENDATION: burn >= 1000
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMMENDATION: samples >= 1000
  samples <- 1000

  # Set the subject names
  names_subjects <- paste(1:n)

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = list(),
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "CONSTANT",
              .subject_names = names_subjects,
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)

  # View the posterior distribution of the number of clusters
  tip_plots$histogram_posterior_number_of_clusters

  # View the trace plot with respect to the posterior number of clusters
  tip_plots$trace_plot_posterior_number_of_clusters

  # Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
  cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

  # If the true labels are available, then show the cluster result via a contigency table
  table(data.frame(true_label = true_labels,
                   cluster_assignment = cluster_assignments))

  # Create the one component graph with minimum entropy
  partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                               .num_components = 1,
                                               .step_size = 0.001)

  # Associate class labels and colors for the plot
  class_palette_colors <- c("Cluster 1" = "blue",
                            "Cluster 2" = 'green',
                            "Cluster 3" = "red")

  # Associate class labels and shapes for the plot
  class_palette_shapes <- c("Cluster 1" = 19,
                            "Cluster 2" = 18,
                            "Cluster 3" = 17)

  # Visualize the posterior similarity matrix by constructing a graph plot of
  # the one-cluster graph. The true labels are used here (below they are not).
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = NA,
                      .subject_class_names = true_labels,
                      .class_colors = class_palette_colors,
                      .class_shapes = class_palette_shapes,
                      .node_size = 2,
                      .add_node_labels = FALSE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  # If true labels are not available, then construct a network plot
  # of the one-cluster graph without any class labels. Also, suppress
  # the subject labels.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = FALSE)
```

