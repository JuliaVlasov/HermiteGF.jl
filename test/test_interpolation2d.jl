@testset "2D" begin
    for i in nvec
        time = @elapsed x = interpolate_2D(:f_3, 
           [:Chebyshev, :Chebyshev], [ep ep], [i i], [Ne Ne], :Hermite, gamma)
           push!(errors["2D"], x[2]) 
           push!(times["2D"], time)
	@test x[2] < max(1.0e-14, 10.0^(-i÷2+1))
    end
end
