# Import the tip library
library(tip)

# Generate a vector of positive integers
# Example: the posterior number of clusters computed after posterior
# sampling in each sampling iteration of the Gibbs sampler.
num_clusters <- c(1,2,2,2,2,3,3,1,2,3,3,3,1,3)

# Generate the plot of the posterior number of clusters versus the
# sampling iteration number in the Gibbs sampler.
ggplot_number_of_clusters_trace(.posterior_number_of_clusters = num_clusters)
