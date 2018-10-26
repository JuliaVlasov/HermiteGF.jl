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
for key in ["1D", "2D", "3D", "4D", "5D"]
   errors[key] = Float64[]
end

@testset "1D" begin
    for i in nvec
        x = interpolate_1D(:f_3, 
        [:Chebyshev], [ep], [i], [Ne], :Hermite, gamma)
        push!(errors["1D"], x[2]) 
	@test x[2] < max(1.0e-14, 10.0^(-i÷2+1))
    end
end
