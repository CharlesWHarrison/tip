# Table Invitation Prior in R (tip)
This R library provides a Gibbs sampler for Bayesian clustering models that utilize the Table Invitation Prior (Harrison, He, and Huang 2022). TIP utilizes pairwise distance information between each observation (i.e. subject) and TIP may be used for clustering vectors, matrices, and higher-order tensors. 

Although the prior used is TIP, there are different options with respect to the likelihood functions. Two options include "CONSTANT" and "NIW". 

The "CONSTANT" likelihood function refers to the scenario where the likleihood function returns a constant value regardless of the input so that likelihood function has no role in the clustering. The "CONSTANT" likelihood may be applied to vector-variate datasets (e.g. the Iris dataset, US Arrests dataset, etc.), matrix-variate datasets (e.g. data pertaining to electroencephalograms (EEGs), grayscale images, etc.), and tensor-variate datasets (i.e. videos, colored-pictures, etc.). 

The "NIW" option refers to the scenario where the likelihood function uses a "Normal-Inverse-Wishart" likelihood function with respect to the current clusters in a given iteration in the Gibbs sampler; the "NIW" option applies only to vectors (i.e. Iris dataset, US Arrests, etc.). 
