import HermiteGF: interpolate_1D
import HermiteGF: interpolate_2D
import HermiteGF: interpolate_3D
import HermiteGF: interpolate_4D
import HermiteGF: interpolate_5D

using Plots
pyplot()

using Test
include("test_interpolation1d.jl")

#=
plot!( nvec, errors["1D"];
       title  = "L2 error scaling",
       markershape = :circle, 
       label  = "1D",
       yscale = :log10)

xlabel!("N")
ylabel!("L2 error")
savefig("errors.png")
=#
