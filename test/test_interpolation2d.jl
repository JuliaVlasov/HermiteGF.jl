f(x, y) = cos.(x.^2 .+ transpose(y.^2))

xmin, xmax, nx = -1, 1, 16
ymin, ymax, ny = -1, 1, 16

# Evaluation points 
nxe, nye = 64, 128
xe = collect(range(xmin, stop=xmax, length=nxe))
ye = collect(range(ymin, stop=ymax, length=nye))
fe = f(xe, ye)
dx = abs(xe[2] - xe[1])
dy = abs(ye[2] - ye[1])

ϵ = 0.1
γ = 3.0

hermite_x = Hermite(Chebyshev(xmax, xmin, nx), ϵ, γ)
hermite_y = Hermite(Chebyshev(ymax, ymin, ny), ϵ, γ)

xk = hermite_x.nodes
yk = hermite_y.nodes
fk = f(xk, yk)

@testset "Hermite 2D" begin
	
    s = interpolate( hermite_x, hermite_y, fk, xe, ye )
    
    max_error = maximum(abs.(s .- f(xe, ye)))
    l2_error = sqrt(trapz((s .- f(xe, ye)).^2, dx, dy))

    println(" Max error $max_error ")
    println(" L2  error $l2_error  ")

    @test l2_error < 1.0e-7

end

radial_x  = Radial(Chebyshev(xmax, xmin, nx), ϵ)
radial_y  = Radial(Chebyshev(ymax, ymin, ny), ϵ)

xk = radial_x.nodes
yk = radial_y.nodes
fk = f(xk, yk)

@testset "Radial 2D" begin

    s = interpolate( radial_x, radial_y, fk, xe, ye )
    
    max_error = maximum(abs.(s .- f(xe, ye)))
    l2_error = sqrt(trapz((s .- f(xe, ye)).^2, dx, dy))

    println(" Max error $max_error ")
    println(" L2  error $l2_error  ")

    @test l2_error < 1.0e-2

end
