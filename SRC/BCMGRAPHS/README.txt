# BCMatch4Graphs

The BCMGRAPHS directory is intended for CLUSTERING and AUGMENTATION of graphs with metadata. 

It includes the following SOURCE FILES:
-bcm_metadata.c : reading the attributes/metadata;
-bcm_augmentation.c : expanding (or augmenting) the original graph by using metadata (algorithm for augmentation by H. Cheng, Y. Zhou, and J.X. Yu., ACM Trans. Knowl. Discov., 2:12.1--12.33, 2011);
-bcm_kmeans.c : clustering of a graph with a standard k-means algorithm;
-bcm_kmeansblock.c : clustering of augmented graphs with attributes with the modified version of the k-means algorithm described in P. D'Ambra, P.S. Vassilevski and L. Cutillo, Applied Mathematics and Computation, Vol. 447, 2023.
-bcm_modularity.c : computing the modularity corresponding to a clustering;
-bcm_entropy.c : computing the entropy wrt the labels corresponding to a clustering of an attributed graph.
 
and the following HEADERS:
bcm_metadata_exp.h : data structures and API for metadata and augmentation.
bcm_graphs.h : data structures and API for bcm_kmeans, bcm_kmeansblock, bcm_modularity and bcm_entropy.
 
