# BCMatch4Graphs rel 1.0

BCMatch4Graphs code is an extension of BootCMatch including functionalities for graph clustering.

It relies on BootCMatch for all the functionalities to embed the graphs in the smooth vectors space related to the graph Laplacian, and to apply a K-means spatial clustering
as described in:  

P. D'Ambra, L. Cutillo, P. S. Vassilevski, Bootstrap AMG for Spectral Clustering. Computational and Mathematical Methods. Vol. 1, 2019, e1020. https://doi.org/10.1002/cmm4.1020.

It also includes the subroutines needed to extend the clustering method to attributed graphs, as described in:

P. D’Ambra, P. S. Vassilevski, L. Cutillo, Extending bootstrap AMG for clustering of attributed graphs, Applied Mathematics and Computation, Vol. 447, 2023. https://doi.org/10.1016/j.amc.2023.127904

All the new functionalities for dealing with graphs are implemented in the new folder SRC/BCMGRAPHS.

See the Test folder for two main programs reading the adjacency matrix of the graph and applying the above methods.
See the Test/RUNS folder for running the tests on sample graphs.

# Installation
```
make clean
make
cd Test
make clean
make

#BCMatch4Graphs
