globalVariables(c("..count.."))

make_invertible <- function(.matrix, .tolerance = 0.01){
  # --- A function used to make a matrix invertible by
  # adding a small number to the matrix's diagonal.
  # If the matrix is invertible, then return the matrix. ---

  # Compute the number of columns
  .dimension <- dim(.matrix)[2]
  if(any(eigen(.matrix)$values <= 0.0)){
    .temp <- .matrix + diag(.dimension)*.tolerance
    if(all(eigen(.matrix)$values > .tolerance)){return(.temp)}
    while(any(eigen(.matrix)$values < 0.0)){
      .temp = .matrix + diag(.dimension)*.tolerance
      .tolerance = .tolerance + .tolerance
    }
    .matrix <- .temp
  }
  return(.matrix)
  # make_invertible(.matrix = matrix(data = 1, nrow = 3, ncol = 3),
  #                 .tolerance = 0.01)
}

#' @export
partition_undirected_graph <-function(.graph_matrix, .num_components, .step_size){

  # --- A function to partition an undirected graph_matrix (i.e. an n x n symmetric matrix) into
  # .num_components components ---

  # 1) .graph_matrix: an n x n matrix where the (i,j)th element denotes
  # the edge weight from vertex i to vertex j

  # 2) .num_components: the desired number of graph components that the
  # graph will be partitioned into

  # 3) .step_size: the transformation max(.graph_matrix[i,j] - .cutoff, 0) is applied
  # iteratively to each element in .graph_matrix until .num_components are obtained
  # (i.e., the components are "islands"). The argument .step_size increments the
  # .cutoff in a loop.

  .flag = 0
  # The initial .cutoff value is set to be very close to the
  # minimum edge weight to speed up the process
  .cutoff = min(.graph_matrix) - .step_size
  while(.flag > -1){

    # Apply the cutoff
    .graph_matrix_temp <- ifelse(.graph_matrix >= .cutoff, .graph_matrix, 0)

    # Partition the graph
    .net <- igraph::graph.adjacency(.graph_matrix_temp,
                                    mode = 'undirected',
                                    weighted = TRUE,
                                    diag = FALSE)

    # Extract the graph components (i.e. clusters)
    .graph_component_members <- igraph::components(.net)$membership


    # Special case when the entire graph with <n> vertices is partitioned into <n> components
    # and the number of desired components is also <n>
    if(length(unique(.graph_component_members)) == dim(.graph_matrix)[1] & .num_components == dim(.graph_matrix)[1]){
      .graph_matrix_temp <- ifelse(.graph_matrix >= .cutoff, .graph_matrix, 0)
      .net <- igraph::graph.adjacency(.graph_matrix_temp,
                                      mode = 'undirected',
                                      weighted = TRUE,
                                      diag = FALSE)
      .graph_component_members <- igraph::components(.net)$membership
      return(list(graph_component_members = .graph_component_members, .cutoff = .cutoff, partitioned_graph_matrix = .graph_matrix_temp))
    }

    # We want the maximum .cutoff value that produces a graph with
    # .num_components, so find the point where the number of graph components
    # is GREATER than the desired number of graph components. Later the
    # cutoff will be decremented to give graph with exactly .num_components.
    # This procedure produces a graph with the desired number of graph components
    # while removing a larger number of (unnecessary) graph edges.

    if(length(unique(.graph_component_members)) > .num_components){
      # Set .cutoff back to optimal .cutoff
      .cutoff <- .cutoff - .step_size

      # Apply the transformation max(0, .graph_matrix[i,j] - .cutoff)
      .graph_matrix_temp <- ifelse(.graph_matrix >= .cutoff, .graph_matrix, 0)

      # Convert the matrix into an igraph graph
      .net <- igraph::graph.adjacency(adjmatrix = .graph_matrix_temp,
                                      mode = 'undirected',
                                      weighted = TRUE,
                                      diag = FALSE)

      # igraph::components(.net)$membership: each subject's island membership
      # cutoff: cutoff s.t. max(.graph_matrix - cutoff,0) gives a graph with .num_components
      # partitioned_graph_matrix: the matrix corresponding to the graph that gives .num_components islands
      return(list(.graph_component_members = igraph::components(.net)$membership,
                  cutoff = .cutoff,
                  partitioned_graph_matrix = .graph_matrix_temp))
    }else{
      # Increase the .cutoff by .step_size
      .cutoff = .cutoff + .step_size
    }

    if(.cutoff >= max(.graph_matrix)){
      # Decrement the .cutoff value by .step_size
      .cutoff <- .cutoff - .step_size

      # Apply the transformation max(0, .graph_matrix[i,j] - .cutoff)
      .graph_matrix_temp <- ifelse(.graph_matrix > .cutoff, .graph_matrix, 0)

      # Convert the matrix into an igraph graph
      .net <- igraph::graph.adjacency(adjmatrix = .graph_matrix_temp,
                                      mode = 'undirected',
                                      weighted = TRUE,
                                      diag = FALSE)

      # igraph::components(.net)$membership: each subject's island membership
      # cutoff: cutoff s.t. max(.graph_matrix - cutoff,0) gives a graph with .num_components
      # partitioned_graph_matrix: the matrix corresponding to the graph that gives .num_components islands

      .graph_component_members <-
        return(list(graph_component_members = igraph::components(.net)$membership,
                    cutoff = .cutoff,
                    partitioned_graph_matrix = .graph_matrix_temp))
    }
  }
  # set.seed(007)
  # partition_graph(.graph_matrix = matrix(abs(rnorm(100)),nrow = 10, ncol = 10),
  #                 .num_components = 3,
  #                 .step_size = 0.01)
}

recode <- function(.posterior_assignments){
  # --- A function to recode the current cluster assignments so that each
  # cluster assignment is in the set {1, 2, 3, ..., K} and there are no gaps.
  # Example: 2, 3, 5, 6, 10, 2, 2, 2, 5 needs to be recoded to 1, 2, 3, 4, 5, 1, 1, 1, 3
  # Example: 1, 2, 3, 4, 5, 2, 2, 2, 5 is recoded to itself since it is already correct ---

  # Compute the current sorted unique cluster assignment values
  .sorted_unique_cluster_values <- sort(unique(.posterior_assignments), decreasing = FALSE)

  # For each current cluster assignment value
  for(j in 1:length(.posterior_assignments)){
    # Replace the current cluster assignment value with the index corresponding to the
    # sorted unique cluster value that the current cluster assignment is equal to
    .posterior_assignments[j] <- which(.sorted_unique_cluster_values == .posterior_assignments[j])
  }
  return(.posterior_assignments)
  # recode(.posterior_assignments = c(2, 3, 5, 6, 10, 2, 2, 2, 5))
}

get_proximity_matrix <- function(.assignments){
  # --- A function to construct a proximity matrix based on
  # a vector of .assignments ---
  # matrix_{i,j} = 1 if subject_i and subject_j belong to the same cluster
  # matrix_{i,j} = 0 if subject_i and subject_j do not belong to the same cluster
  return(outer(.assignments, .assignments, function(x, y) as.integer(x==y)))
}
