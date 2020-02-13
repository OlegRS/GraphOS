# GraphOS
C++ library for network theory simulations

# Purpose
The library emerged from my studies in various aspects of network theory. The internal representation of the graphs allows working with large sparse networks in a memory efficient way. Can be useful as a tool for working with random graphs and various spin models on networks, especially if the networks are sparse.

# Functionality
## Generation
GraphOS can be used to generate random graphs from various user defined ensembles by means of Metropolis dynamics.

## Analysis
Some local and global features of graphs (such as degrees, clustering coefficients, modularity, etc.) can be computed using GraphOS.

## Spin models on networks
Various spin models (such as Spin Glasses, Boltzmann Machines, Multi-component Curie-Weiss model, etc.) can be simulated in GraphOS.

## Classification
GraphOS can be used to play with unsupervised/semisupervised classification of nodes. Currently only the brute-force precise classifiers are implemented (one based on modularity and another on likelihood maximisation), so classification is only feasible for very small (~10 nodes) graphs.

## Communication with other software
It is possible to save the graphs in .sif-.attrs format recognised by cytoscape (https://cytoscape.org).

# Usage
Use cases can be found in separate repositories dedicated to particular research projects on network theory and random matrices.

# Things to improve
A lot of things can be added and improved in GraphOS. The code should at least be refactored and having a Python interface would also be useful.
