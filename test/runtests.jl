using HermiteGF

using Test

include("trapz.jl")
include("ndgrid.jl")
include("test_interpolation1d.jl")
include("test_interpolation2d.jl")

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
