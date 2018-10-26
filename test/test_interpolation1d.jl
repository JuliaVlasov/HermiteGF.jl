@testset "1D" begin
    for i in nvec
        time = @elapsed x = interpolate_1D(:f_3, 
        [:Chebyshev], [ep], [i], [Ne], :Hermite, gamma)

        push!(times["1D"], time)
        push!(errors["1D"], x[2]) 
	@test x[2] < max(1.0e-14, 10.0^(-i÷2+1))
    end
end
