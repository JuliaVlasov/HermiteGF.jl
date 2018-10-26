###
# In this file the test for the error scaling with the dimension is implemented.
# In particular, the dependence of the L2 error from the number of points per dimension.
###
import HermiteGF: interpolate_1D
import HermiteGF: interpolate_2D

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

@testset "1D" begin
    for i in nvec
        time = @elapsed x = interpolate_1D("f_3", [:Chebyshev], [ep], [i], [Ne], :Hermite, gamma)
        push!(times["1D"], time)
        push!(errors["1D"], x[2]) 
	println(x[1])
	@test x[2] < max(1.0e-14, 10.0^(-i÷2+1))
    end
end


for i in nvec
  time = @elapsed x = interpolate_2D("f_3", ["Chebyshev", "Chebyshev"], [ep ep], [i i], [Ne Ne], "Hermite", gamma);
  push!(errors["2D"], x[2]) 
  push!(times["2D"], time)
end

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

  #=
  errors_5D[:, ind] = interpolate_5D("f_3", ["Chebyshev","Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep ep], [i i i i i], [Ne Ne Ne Ne Ne], "Hermite", gamma);
  errors_4D[:, ind] = interpolate_4D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep], [i i i i], [Ne Ne Ne Ne], "Hermite", gamma);=#
  #errors_3D[:, ind] = interpolate_3D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep], [i i i], [Ne Ne Ne], "Hermite", gamma);
#plot(nvec, errors_4D[1, :], marker="o", label="4D");
plot(nvec, errors_3D[1, :], marker="<", label="3D");
plot(nvec, errors_2D[1, :], marker="<", label="2D");
plot(nvec, errors_1D[1, :], marker=">", label="1D");
ax = gca() # Get the handle of the current axis
ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
legend(loc="upper right",fancybox="true")
xlabel("N")
ylabel("Maximum error")=#

# Plot the L2 error
#plot(nvec, errors_4D[2, :], marker="o", label="4D");
#plot(nvec, errors_3D[2, :], marker="<", label="3D");
#plot(nvec, errors_2D[2, :], marker="<", label="2D");
