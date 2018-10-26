# HermiteGF.jl

[![Build Status](https://travis-ci.org/JuliaVlasov/HermiteGF.jl.svg?branch=master)](https://travis-ci.org/JuliaVlasov/HermiteGF.jl)

Stable Gaussian radial basis function interpolation based on HermiteGF expansion

Author: Anna Yurova

This is an implementation of the method described in the paper

[*STABLE EVALUATION OF GAUSSIAN RADIAL BASIS FUNCTIONS USING HERMITE POLYNOMIALS*](https://arxiv.org/abs/1709.02164)

by Anna Yurova and Katharina Kormann.

- HermiteGF-tensor in 1-5D. The implementation of 4-5D cases is parallel. 5D tests have to be run on a cluster to be finished in a reasonable time.

Installation

```julia
using Pkg
Pkg.add("https://gitlab.mpcdf.mpg.de/clapp/HermiteGF.jl")
using HermiteGF
```
