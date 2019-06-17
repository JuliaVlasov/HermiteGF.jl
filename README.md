# HermiteGF.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.org/JuliaVlasov/HermiteGF.jl.svg?branch=master)](https://travis-ci.org/JuliaVlasov/HermiteGF.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/HermiteGF.jl/latest)
[![codecov](https://codecov.io/gh/JuliaVlasov/HermiteGF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/HermiteGF.jl)

Stable Gaussian radial basis function interpolation based on HermiteGF expansion

This is an implementation of the method described in the paper

Anna Yurova and Katharina Kormann [*Stable evaluation of guassian radial basis functions using Hermite polynomials*](https://arxiv.org/abs/1709.02164).

Installation

```julia
using Pkg
Pkg.add("https://gitlab.mpcdf.mpg.de/clapp/HermiteGF.jl")
using HermiteGF
```
