# --- Functions used to generate useful plots for Bayesian clustering ---
# library(ggplot2)
# library(igraph)
# library(network)
# library(GGally)
# library(network)

#' @title Plot connected points using ggplot2.
#' @description A function to that produces a ggplot2 plot of .y versus .x
#' where points are added via geom_point() and the points are connected via geom_line().
#' @param .x The variable on the horizontal axis.
#' @param .y The variable on the vertical axis.
#' @param .ylab The label on the vertical axis.
#' @param .xlab The label on the horizontal axis.
#' @importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab
#' @importFrom rlang .data
#' @export
ggplot_line_point <- function(.x, .y, .xlab = "", .ylab = ""){
  # --- A function to plot a line and the corresponding points ---
  plot_df <- data.frame(x = .x, y = .y)
  p <- ggplot(data = plot_df, aes(x = .data$x, y = .data$y)) + geom_line()
  p <- p + geom_point() + xlab(.xlab) + ylab(.ylab)
  return(p)
}

#' @title Plot the posterior distribution of the number of clusters.
#' @description A function that produces a ggplot bar chart (i.e. geom_bar) that corresponds
#' to the posterior number of clusters. The vertical axis is normalized so that it displays
#' the posterior probability.
#' @param .posterior_number_of_clusters A vector of the number of clusters after posterior sampling
#' for each iteration in the Gibbs sampler.
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab scale_x_continuous
#' @importFrom rlang .data
#' @export
ggplot_number_of_clusters_hist <- function(.posterior_number_of_clusters){
  .table = table(.posterior_number_of_clusters)
  # --- A function to construct the posterior distribution of the number of clusters ---
  .df_tip <- data.frame(y = as.numeric(.table)/sum(.table),
                        x = as.numeric(names(.table)))
  plot <- ggplot(data = .df_tip) + geom_bar(aes(x = .data$x, y = .data$y), stat = "identity")
  plot <- plot + xlab("Number of Clusters")
  plot <- plot + ylab("Posterior Probability")
  plot <- plot + scale_x_continuous(breaks = 1:max(.df_tip$x))
  return(plot)
}

#' @title Plot the trace plot of the posterior number of clusters.
#' @description A function that produces a ggplot2 trace plot (i.e. geom_line)
#' with respect to the posterior number of clusters.
#' @param .posterior_number_of_clusters A vector of the number of clusters after posterior sampling
#' for each iteration in the Gibbs sampler.
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab
#' @importFrom rlang .data
#' @export
ggplot_number_of_clusters_trace <- function(.posterior_number_of_clusters){
  # --- A function to construct the posterior distribution of the number of clusters ---
  .df_tip <- data.frame(iteration = 1:length(.posterior_number_of_clusters),
                        num_clusters = .posterior_number_of_clusters)
  plot <- ggplot(data = .df_tip, aes(y = .data$num_clusters, x = .data$iteration)) + geom_line()
  plot <- plot + xlab("Gibbs Sampling Iteration")
  plot <- plot + ylab("Posterior Number of Clusters")
  return(plot)
}

#' @title Visualize the posterior similarity matrix (i.e. posterior probability matrix)
#' @description A function that produces a ggnet2 network plot to visualize the posterior similarity matrix (i.e. the matrix of posterior probabilities).
#' @param .matrix_graph A matrix M where each element Mij corresponds to the posterior
#' probability that subjects i and j are in the same cluster.
#' @param .subject_names An optional vector of subject names that will appear in the graph plot.
#' @param .subject_class_names An optional vector of class names that will influence each vertex's color and shape.
#' @param .class_colors An optional named vector of colors. The vector names are required to be
#' the unique .subject_class_names whereas the vector values are required to be the colors.
#' @param .class_shapes An optional named vector of shapes. The vector names are required to be
#' the unique .subject_class_names whereas the vector values areq required to be integers
#' (i.e. pch values like 15, 16, 17, and so on).
#' @param .random_seed The plot uses the random layout, so set a seed for reproducibility.
#' @param .node_size The size of each node (i.e. vertex) in the graph plot.
#' @param .add_node_labels TRUE or FALSE. Should individual node labels be added to each node (i.e. vertex) in the graph plot?
#' @import GGally
#' @import network
#' @export
ggnet2_network_plot <- function(.matrix_graph, .subject_names = vector(), .subject_class_names = NA,
                             .class_colors, .class_shapes, .random_seed = 007, .node_size = 6,
                             .add_node_labels = TRUE){
  # --- A function to construct a network plot ---

  # Construct the network object
  .network_temp <- network::network(x = .matrix_graph,
                                    directed = FALSE,
                                    ignore.eval = FALSE,
                                    names.eval = "weights")

  # Add labels to the graph nodes (i.e., the subjects)
  if(length(.subject_names) != dim(.matrix_graph)[1]){
    # network::network.vertex.names(.network_temp) <- paste("Subject", 1:dim(.matrix_graph)[1], sep = "")
  }else{
    network.vertex.names(.network_temp) <- .subject_names
  }

  if(length(.subject_class_names) == dim(.matrix_graph)[1]){
    .network_temp %v% "Category" <- as.character(.subject_class_names)

    # Set a random seed so that the graph node (subject) layouts are reproducible
    set.seed(.random_seed)

    ggnet2(.network_temp,
           shape = "Category",
           shape.palette = .class_shapes,
           color = "Category",
           color.palette = .class_colors,
           label = .add_node_labels,
           node.size = .node_size)
  }else{
    ggnet2(.network_temp, label = .add_node_labels)
  }
}

