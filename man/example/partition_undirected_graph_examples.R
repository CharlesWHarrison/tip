# Import the tip library
library(tip)

# Choose an arbitrary random seed to generate the data
set.seed(4*8*15*16*23*42)

# Generate a symmetric posterior probability matrix
# Each element is the probability that the two subjects belong
# to the same cluster
n1 <- 10
posterior_prob_matrix <- matrix(NA, nrow = n1, ncol = n1)
for(i in 1:n1){
  for(j in i:n1){
    if(i != j){
      posterior_prob_matrix[i,j] <- runif(n=1,min=0,max=1)
      posterior_prob_matrix[j,i] <- posterior_prob_matrix[i,j]
    }else{
      posterior_prob_matrix[i,j] <- 1.0
    }
  }
}

# Generate a one-cluster graph (i.e., partitioned_graph_matrix)
partition_undirected_graph(.graph_matrix = posterior_prob_matrix,
                           .num_components = 1,
                           .step_size = 0.001)
