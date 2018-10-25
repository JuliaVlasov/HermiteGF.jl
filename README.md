# HermiteGF.jl

Stable Gaussian radial basis function interpolation based on HermiteGF expansion

HermiteGF stabilization code for the RBF interpolation. 

Author: Anna Yurova

This is an implementation of the method described in the paper

"STABLE EVALUATION OF GAUSSIAN RADIAL BASIS FUNCTIONS USING HERMITE POLYNOMIALS"

by Anna Yurova and Katharina Kormann.

https://arxiv.org/abs/1709.02164

- HermiteGF-tensor in 1-5D. The implementation of 4-5D cases is parallel. 5D tests have to be run on a cluster to be finished in a reasonable time.

In order to install Julia on your computer, perform the following steps:

- Download julia from https://julialang.org/downloads/ and unzip it.  Note: The plotting currently only works with Julia v0.5.
- Note that you need to make sure curl and cmake are installed. On Ubuntu:
  sudo apt-get install curl
  sudo apt-get install cmake
- Download atom from https://atom.io/ and install it (On Ubuntu the package manager can be used).
- Install uber-juno through installation manager in atom.
- Set the Julia path to the Julia binary that was installed in the first step (Use settings -> packages -> Julia -> setting).
- Start Julia.

Before running simulations. It is necessary to install the package

```julia
using Pkg
Pkg.clone("https://gitlab.mpcdf.mpg.de/clapp/HermiteGF.jl")
using HermiteGF
```

In order to run parallel simulations, it is necessary to start Julia with appropriate amount of processes. Command line example:

```
julia -p 16 -L $HOME/hermiteGF/Julia/init.jl $HOME/hermiteGF/Julia/test_dependence_on_N.jl
```

We ask you to cite the following reference in scientific publica-
tions which contain results obtained with this software and developments:
*A. Yurova, K. Kormann
“Stable evaluation of Gaussian radial basis functions using Hermite polyno-
mials”*
