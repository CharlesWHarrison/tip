# Clustering vectors, matrices, and tensors using the Table Invitation Prior (TIP) in R 
``` 
# Install from CRAN
install.packages("tip") 
```

```
# Install from GitHub
devtools::install_github("STATS-ML/tip")
```

This R library provides a Gibbs sampler for Bayesian clustering models that utilize the Table Invitation Prior (TIP) introduced by Harrison, He, and Huang (2022). TIP utilizes pairwise distance and pairwise similarity information between the observed data (i.e. subjects). The term ''subject'' is used to refer to an individual vector, matrix, or higher-order tensors. 
1. **Vector-variate subject example**: in the Iris dataset there are 150 observed flowers, and each flower's characteristics (not including their species) is captured by a ``4 x 1`` vector. Each flower is considered as an individual subject, so there are 150 subjects.
2. **Matrix-variate subject example**: a single X-ray is taken for 57 adults, and each X-ray image is stored as a ``512 x 512`` matrix where each value in each matrix varies between zero and one (i.e., a grayscale image). Each X-ray is considered as an individual subject, so there are 57 subjects.
3. **Tensor-variate subject example**: 23 adults have an fMRI taken. Each fMRI image corresponds to a 3-way tensor, each of the 23 3-way tensors correspond to an individual subject, so there are 23 subjects.

In Bayesian clustering, the goal is sample from the posterior distribution which is given by 
$$P(\mathbf{c}|\mathbf{x}) \propto P(\mathbf{X}|\mathbf{c})P(\mathbf{c})$$ 
where $P(\mathbf{X}|\mathbf{c})$ is the likelihood function and $P(\mathbf{c})$ is the prior distribution. In this case, $P(\mathbf{c})$ refers to the Table Invitation Prior (TIP). TIP is flexible and may be incorporated with a variety of different likelihood functions designed for different types of data. The current options for the likliehood function are the following:

1. The ``.likelihood_model = "CONSTANT"`` is the fastest option and can be used for vectors, matrices, and higher-order tensors (i.e., ``.data`` is not used). The "CONSTANT" option returns a constant likelihood function value regardless of the observed data so that likelihood function has no role in the clustering. The "CONSTANT" likelihood option may be used for vector-variate datasets (e.g. the Iris dataset, US Arrests dataset, etc.), matrix-variate datasets (e.g. data pertaining to electroencephalograms (EEGs), grayscale images, etc.), and higher-order tensor-variate datasets (i.e. videos, colored-pictures, etc.). 

2. The ``.likelihood_model = "NIW"`` option can be used for vector-datasets only (i.e., ``.data`` is a ``data.frame``). The "NIW" option uses a "Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of vector-variate data include the Iris dataset, US Arrests dataset, and so on.  

3. The ``.likelihood_model = "MNIW"`` option can be used for matrix-variate data (i.e., ``.data`` is a ``list`` of matrices). The "MNIW" option uses a "Matrix-Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of matrix-variate data include EEG data, grayscale images (i.e., X-rays or black and white photographs), graph data, and so on. Note that "MNIW" may be used for vector-variate data by passing a list of matrices that have either 1 row and multiple columns or matrices that have one column and multiple rows.

## Clustering the Iris Dataset (i.e. vectors) with a Normal-Inverse-Wishart (NIW) likelihood and a TIP prior
```
##### BEGIN EXAMPLE 1: Clustering the Iris Dataset (vector clustering) #####
  ##### Prior Distribution: Table Invitation Prior (TIP)
  ##### Likelihood Model: Normal Inverse Wishart (NIW)

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
  # RECOMMENDATION: burn >= 1000
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
  # Note: Subject labels may be suppressed using .add_node_labels = FALSE.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  ##### END EXAMPLE 1: Vector Clustering (NIW) #####



  ##### BEGIN EXAMPLE 2: Clustering the US Arrests Dataset (vector clustering) #####
  ##### Prior Distribution: Table Invitation Prior (TIP)
  ##### Likelihood Model: Normal Inverse Wishart (NIW)

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
  # Note: Subject labels may be suppressed using .add_node_labels = FALSE.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  ##### END EXAMPLE 2: Clustering the US Arrests Dataset (vector clustering) #####



  ##### BEGIN EXAMPLE 3: Matrix Clustering (MNIW) #####
  ##### Prior Distribution: Table Invitation Prior
  ##### Likelihood Model: Matrix Normal Inverse Wishart (MNIW)

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

  # Initialize a matrix of NA values
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
  names_subjects <- paste(1:dim(iris)[1])

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
  # Note: Subject labels may be suppressed using .add_node_labels = FALSE.
  ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                      .subject_names = names_subjects,
                      .node_size = 2,
                      .add_node_labels = TRUE)

  ##### END EXAMPLE 3: Matrix Clustering (MNIW) #####
```
