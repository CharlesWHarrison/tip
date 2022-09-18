# # Import the TIP source code
# # source("tip.R")
# # source("plot_fn.R")
#
# # Import the iris dataset
# data(iris)
#
# # The first 4 columns are the data whereas
# # the 5th column refers to the true labels
# X <- data.matrix(iris[,c("Sepal.Length",
#
#                          "Sepal.Width",
#                          "Petal.Length",
#                          "Petal.Width")])
#
# # Extract the true labels
# true_labels <- iris[,"Species"]
#
# # Compute the similarity matrix
# distance_matrix <- data.matrix(dist(iris[,1:4]))
# temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])
#
# # Compute the point estimate for the nearest neighbors using
# # univariate multiple change point detection (i.e.)
# init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)
#
# # Set the number of burn-in iterations in the Gibbs samlper
# burn <- 100
#
# # Set the number of sampling iterations in the Gibbs sampler
# samples <- 100
#
# # Run TIP clustering with a Normal-Inverse-Wishart likelihood function
# tip1 <- tip(.data = data.matrix(iris[,1:4]),
#             .burn = burn,
#             .samples = samples,
#             .similarity_matrix = exp(-1.0*temperature*distance_matrix),
#             .init_num_neighbors = init_num_neighbors,
#             .likelihood_model = "NONE",
#             .subject_names = paste("Flower",1:dim(iris)[1]),
#             .num_cores = 1)
#
# # Show the posterior distribution of the number of clusters
# ggplot_number_of_clusters_hist(.posterior_number_of_clusters = tip1$posterior_number_of_clusters)
#
# # Show the trace plot of the posterior number of clusters
# ggplot_number_of_clusters_trace(.posterior_number_of_clusters = tip1$posterior_number_of_clusters)
#
# # Extract posterior cluster assignments using PEAR
# cluster_assignments <- mcclust::maxpear(psm = tip1$posterior_similarity_matrix)$cl
#
# # Create the one component graph with minimum entropy
# partition_list <- partition_undirected_graph(.graph_matrix = tip1$posterior_similarity_matrix,
#                                              .num_components = 1,
#                                              .step_size = 0.001)
#
# # Associate class labels and colors for the plot
# class_palette_colors <- c("setosa" = "blue",
#                           "versicolor" = 'green',
#                           "virginica" = "orange")
#
# # Associate class labels and shapes for the plot
# class_palette_shapes <- c("setosa" = 19,
#                           "versicolor" = 18,
#                           "virginica" = 17)
#
# # Visualize the posterior similarity matrix
# # by constructing the one-cluster graph and with the one-cluster graph using a network plot
# ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
#                 .subject_names = NA,
#                 .subject_class_names = true_labels,
#                 .class_colors = class_palette_colors,
#                 .class_shapes = class_palette_shapes,
#                 .node_size = 2,
#                 .add_node_labels = FALSE)
