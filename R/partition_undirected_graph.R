#' @title Partition an undirected graph
#' @description A function that iteratively applies the transformation max(0, .graph_matrix - cutoff) until
#' there are <.num_components> graph components where cutoff = cutoff + .step_size. This is used to generate the one-cluster graph and plot.
#' @param .graph_matrix Matrix: a symmetric matrix that the analyst wishes to decompose into <.num_components> components.
#' @param .num_components Positive integer: the number of components that the analyst wishes to decompose <.graph_matrix> into.
#' @param .step_size Positive numeric: the size of the update for the cutoff in the transformation max(0, .graph_matrix - cutoff)
#' where cutoff = cutoff + .step_size.
#' @returns List with three elements:
#' 1) graph_component_members: a vector of graph component assignments;
#' 2) cutoff: the value max(0, g_{i,j} - cutoff) so that there are <.num_components> components in the graph
#' 3) partitioned_graph_matrix: the graph with <.num_components> components (parts).
#' @example example/partition_undirected_graph_examples.R
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
