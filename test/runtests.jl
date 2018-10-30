using HermiteGF

using Test

include("trapz.jl")
include("test_interpolation1d.jl")
include("test_interpolation2d.jl")
include("test_interpolation3d.jl")
include("test_interpolation4d.jl")
include("test_interpolation5d.jl")

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
