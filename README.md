# GraphOS
C++ library for network theory simulations

# Purpose
The library emerged from my studies in various aspects of network theory. It is far from being any comprehensive, but rather simple and, thus, allowing the user to be in full control of what it does. The internal representation of graphs allows working with the large sparse networks in a memory efficient way.

# Functionality
## Generation
GraphOS can be used to generate random graphs from various user defined ensembles by means of the MCMC algorithms. For example, each node of the network can carry an attribute (label), and each pair of nodes can be connected with probability that depends on the corresponding pair of node labels. GraphOS was also used to sample from various exponential network ensembles using Metropolis dynamics.

## Analytics
Some local and global features of graphs (such as degree, degree distribution, modularity, etc.) can be computed using GraphOS.

## Spin models on networks
Various spin models (such as Spin Glasses, Boltzmann Machines, Multi-component Curie-Weiss model, etc.) can be simulated in GraphOS.

## Classification
GraphOS can be used to play with unsupervised/semisupervised classification of nodes. Currently only the precise brute-force classifiers are implemented (one based on modularity and another based on likelihood maximisation), so classification is only feasible for small (~10 nodes) graphs.

## Communication with other software
It is possible to save the graphs in .sif-.attrs format recognised by cytoscape (https://cytoscape.org).

# Usage
Use cases can be found in separate repositories dedicated to particular research projects on network theory and random matrices.

# Things to improve
A lot of things can be added and improved in GraphOS.
For example, it would be nice to implement Maximum Likelihood classification of nodes with some good heuristic (at least greedy heuristic) and see if it can compare with spectral clustering or other good algorithms of community detection. 
It would also be nice to make a Python interface for GraphOS.
