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
2. **Matrix-variate subject example**: a single X-ray is taken for 57 adults, and each X-ray image is stored as a ``512 x 512`` matrix where each value in the matrix varies between zero and one (i.e., a grayscale image). Each X-ray is considered as an individual subject, so there are 57 subjects.
3. **Tensor-variate subject example**: 23 adults have an fMRI taken. Each fMRI image corresponds to a 3-way tensor, each of the 23 3-way tensors correspond to an individual subject, so there are 23 subjects.

Although the prior used is TIP, there are different options with respect to the likelihood functions. Currently there are three options for the likelihood model:

1. The ``.likelihood_model = "CONSTANT"`` is the fastest option and can be used for vectors, matrices, and higher-order tensors (i.e., ``.data`` is not used). The "CONSTANT" option returns a constant likelihood function value regardless of the observed data so that likelihood function has no role in the clustering. The "CONSTANT" likelihood option may be used for vector-variate datasets (e.g. the Iris dataset, US Arrests dataset, etc.), matrix-variate datasets (e.g. data pertaining to electroencephalograms (EEGs), grayscale images, etc.), and higher-order tensor-variate datasets (i.e. videos, colored-pictures, etc.). 

2. The ``.likelihood_model = "NIW"`` option can be used for vector-datasets only (i.e., ``.data`` is a ``.data.frame``). The "NIW" option uses a "Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of vector-variate data include the Iris dataset, US Arrests dataset, and so on.  

3. The ``.likelihood_model = "MNIW"`` option can be used for matrix-variate data (i.e., ``.data`` is a ``list`` of matrices). The "MNIW" option uses a "Matrix-Normal-Inverse-Wishart" likelihood function with respect to the current clusters and the new cluster in a given iteration in the Gibbs sampler. Examples of matrix-variate data include EEG data, grayscale images (i.e., X-rays or black and white photographs), graph data, and so on. Note that "MNIW" may be used for vector-variate data by passing a list of matrices that have either 1 row and multiple columns or matrices that have one column and multiple rows.

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
burn <- 1000

# Set the number of sampling iterations in the Gibbs sampler
# A very good result for Iris may be obtained by setting samples <- 1000
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
