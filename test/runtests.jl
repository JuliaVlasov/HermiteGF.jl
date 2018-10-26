import HermiteGF: interpolate_1D
import HermiteGF: interpolate_2D
import HermiteGF: interpolate_3D
import HermiteGF: interpolate_4D
import HermiteGF: interpolate_5D

using Plots
pyplot()

using Test

# Number of collocation point per dimension
nvec  = collect(5:35)
# Value of the shape parameter (in this test the same in all directions, but can be different)
ep    = 0.1
# Number of evaluation points per dimension
Ne    = 53 
# scaling parameter gamma
gamma = 3 

# Initializing error vectors
# first row - maximum error
# second row - L2 error
errors = Dict()
times = Dict()
for key in ["1D", "2D", "3D", "4D", "5D"]
   errors[key] = Float64[]
   times[key]  = Float64[]
end

include("test_interpolation1d.jl")
include("test_interpolation2d.jl")
#include("test_interpolation3d.jl")
#include("test_interpolation4d.jl")
#include("test_interpolation5d.jl")

p = plot(title  = "L2 error scaling")

for dim in keys(errors)
    plot!(p, nvec, errors[dim];
          markershape = :circle, 
	  label  = string(dim),
          yscale = :log10)
end

display(p)
xlabel!("N")
ylabel!("L2 error")
savefig("errors.png")

for (key, value) in sort(collect(times), by=last)
    println(rpad(key, 25, "."), lpad(round(value, digits=1), 6, "."))
end
