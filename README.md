# Clustering vectors, matrices, and tensors using the Table Invitation Prior (TIP) in R 
This R library provides a Gibbs sampler for Bayesian clustering models that utilize the Table Invitation Prior (Harrison, He, and Huang 2022). TIP utilizes pairwise distance information between each observation (i.e. subject) and TIP may be used for clustering vectors, matrices, and higher-order tensors. 

Although the prior used is TIP, there are different options with respect to the likelihood functions. Currently there are two options for the likelihood model:

1. The ```.likelihood_model = "CONSTANT"``` option refers to the scenario where the likleihood function returns a constant value regardless of the input (i.e. vectors, matrices, tensors) so that likelihood function has no role in the clustering. The "CONSTANT" likelihood may be applied to vector-variate datasets (e.g. the Iris dataset, US Arrests dataset, etc.), matrix-variate datasets (e.g. data pertaining to electroencephalograms (EEGs), grayscale images, etc.), and tensor-variate datasets (i.e. videos, colored-pictures, etc.). 

2. The ```.likelihood_model = "NIW"``` option refers to the scenario where the likelihood function uses a "Normal-Inverse-Wishart" likelihood function with respect to the current clusters in a given iteration in the Gibbs sampler. The "NIW" model is applicable to vector-variate datasets only (i.e. the Iris dataset, US Arrests dataset, etc.).  


## Clustering the Iris Dataset (i.e. vectors) with a Normal-Inverse-Wishart (NIW) likelihood and a TIP prior
```
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
# A very good result for Iris may be obtained by setting burn <- 1000
burn <- 50

# Set the number of sampling iterations in the Gibbs sampler
# A very good result for Iris may be obtained by setting samples <- 1000
samples <- 50

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

# Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

# If the true labels are available, then show the cluster result via a contigency table
table(data.frame(true_label = true_labels,
                 cluster_assignment = cluster_assignments))
                 
# Produce plots for the Bayesian Clustering Model
plot(tip1) 

# Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

# If the true labels are available, then show the cluster result via a contigency table
table(data.frame(true_label = true_labels,
                 cluster_assignment = cluster_assignments))
```
