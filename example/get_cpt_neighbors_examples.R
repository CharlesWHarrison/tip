# Import the tip library
library(tip)

# Choose an arbitrary random seed
set.seed(007)

# Generate some data (i.e., 20 subjects described by a 5 x 1 vector)
X <- matrix(rnorm(10*10),nrow=20,ncol=5)

# Compute the pairwise distances between the subjects
distance_matrix <- data.matrix(dist(X))

# For each subject, find the estimate for the number of similar subjects
get_cpt_neighbors(.distance_matrix = distance_matrix)
