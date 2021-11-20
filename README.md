# MOMPKM
Codes for Median of Means Power k-means. Detailed paper, "Uniform Concentration Bounds toward a Unified  Framework for Robust Clustering"  to appear in NeurIPS'21

All the functions are in 'functions.R'. An example run is given in the pdf and RMD files.

Description of MOMPKM:

Implements MOMPKM described in the paper.

# Inputs: 

X 	: n * p data matrix whose rows denote the observations.

k 	: No. of clusters.

L 	: No. of partitions.

s 	: Initial value for parameter s in power k-means. Default value is -1.

eta     : Rate of increase of s (s= s * eta). Default value is 1.02.

alpha   : Learning rate for Adagrad Default is 0.1.

verbose : If TRUE, prints iteration numbers. Default value is FALSE.

tmax    : Maximum number of iterations to run the algorithm. Default is 100.

# Outputs:

label   : Class label of the n data points, returned as an n length vector.

theta   : A k * p matrix, whose rows denote the cluster centroids.   
