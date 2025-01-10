# Elliptic_Problems_Series

MATLAB code for nonparametric Bayesian inference in elliptic PDEs with Gaussian series priors.

Author: Matteo Giordano, https://matteogiordano.weebly.com.

This repository is associated with the article "Bayesian elliptic inverse problems with Gaussian series priors" by Matteo Giordano. The abstract of the paper is:

"This work considers the inverse problem of recovering the unknown conductivity function in an elliptic partial differential equation (PDE) from interior measurements of the PDE solution. Nonparametric Bayesian procedures with Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis are devised. Using recent results from the theory of Bayesian inverse problems, the associated posterior distributions are shown to be asymptotically consistent as the sample size increases, leading to posterior mean estimates with asymptotically vanishing reconstruction errors. An implementation of posterior-based inference is provided via a robust Markov chain Monte Carlo method for functions, and illustrated in a numerical
simulation study."

This repository contains the MATLAB code to replicate the numerical simulation study presented in Section 3 of the article. It contains the following:

GenerateObservations.m code to generate the observations (discrete point evaluations of an elliptic PDE solution corrupted by additive Gaussian measurement errors)
pCNSeries.m code to implement posterior inference based on Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis, via the pCN algorithm.

For questions or for reporting bugs, please e-mail Matteo Giordano (matteo.giordano@unito.it)

Please cite the following article if you use this repository in your research: Giordano, M (2025). Bayesian elliptic inverse problems with Gaussian series priors
