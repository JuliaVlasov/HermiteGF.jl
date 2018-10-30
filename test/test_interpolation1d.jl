using Printf

@testset "Hermite 1D" begin

    for nx in 5:5:35

        ϵ, γ = 0.1, 3.0
        xmin, xmax = -1, 1

        interp = Hermite(Chebyshev(xmin, xmax, nx), 0.1, 3.0)    

        xk = interp.nodes
        fk = cos.(xk.^2)

	nxe = 128
        xe  = collect(range(xmin, stop=xmax, length=nxe))
        dx  = xe[2]-xe[1]
        fe  = cos.(xe.^2)

        s   = interpolate( interp, fk, xe )

        l2_error = sqrt(trapz((s - fe).^2, dx))
        l1_error = maximum(abs.(s .- fe))
	@printf( "%02d points - L2 : %.3e - L∞ : %.3e \n",
	        nx, l2_error, l1_error )

	@test l2_error < max(1.0e-14, 10.0^(-nx÷2+1))
    end

end

@testset "Radial 1D" begin

    for nx in 5:5:35

        ϵ, γ = 0.1, 3.0
        xmin, xmax = -1, 1

        interp = Radial(Chebyshev(xmin, xmax, nx), 0.1)    

        xk = interp.nodes
        fk = cos.(xk.^2)

	nxe = 128
        xe  = collect(range(xmin, stop=xmax, length=nxe))
        dx  = xe[2]-xe[1]
        fe  = cos.(xe.^2)

        s   = interpolate( interp, fk, xe )

        l2_error = sqrt(trapz((s - fe).^2, dx))
        l1_error = maximum(abs.(s .- fe))
	@printf( "%02d points - L2 : %.3e - L∞ : %.3e \n",
	        nx, l2_error, l1_error )

	#@test l2_error < max(1.0e-14, 10.0^(-nx÷2+1))
    end

end
