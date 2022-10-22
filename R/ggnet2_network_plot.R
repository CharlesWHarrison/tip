#' @title Visualize the posterior similarity matrix (i.e., posterior probability matrix)
#' @description A function that produces a ggnet2 network plot to visualize the posterior similarity matrix (i.e., the matrix of posterior probabilities).
#' @param .matrix_graph Matrix. A matrix M where each element Mij corresponds to the posterior
#' probability that subjects i and j are in the same cluster.
#' @param .subject_names Vector of characters. An optional vector of subject names that will appear in the graph plot.
#' @param .subject_class_names A vector of characters. An optional vector of class names that will influence each vertex's color and shape.
#' @param .class_colors Named vector of characters. An optional named vector of colors. The vector names are required to be
#' the unique .subject_class_names whereas the vector values are required to be the colors.
#' @param .class_shapes Named vector of integers. An optional named vector of shapes. The vector names are required to be
#' the unique .subject_class_names whereas the vector values are required to be positive integers
#' (i.e., pch values like 15, 16, 17, and so on).
#' @param .random_seed Numeric. The plot uses the random layout, so set a seed for reproducibility.
#' @param .node_size Positive integer. The size of each node (i.e., vertex) in the graph plot.
#' @param .add_node_labels Boolean. TRUE or FALSE. Should individual node labels be added to each node (i.e., vertex) in the graph plot?
#' @returns ggnet2 plot. A network plot with respect to the undirected network given by .matrix_graph. This is used to visualize the posterior similarity matrix.
#' @import GGally
#' @import network
#' @export
ggnet2_network_plot <- function(.matrix_graph, .subject_names = vector(), .subject_class_names = vector(),
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
