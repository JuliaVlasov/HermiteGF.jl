@testset "Hermite 3D" begin

    f(x, y, z) = cos(x^2 + y^2 + z^2)

    xmin, xmax, nx = -1, 1, 32
    ymin, ymax, ny = -1, 1, 32
    zmin, zmax, nz = -1, 1, 32

    ϵ = 0.1
    γ = 3.0

    hermite_x = Hermite(Chebyshev( xmin, xmax, nx ), ϵ, γ)
    hermite_y = Hermite(Chebyshev( ymin, ymax, ny ), ϵ, γ)
    hermite_z = Hermite(Chebyshev( zmin, zmax, nz ), ϵ, γ)

    xk = hermite_x.nodes
    yk = hermite_y.nodes
    zk = hermite_z.nodes

    nxe, nye, nze = 32, 64, 128

    xe = collect(range(xmin, stop=xmax, length=nxe))
    ye = collect(range(ymin, stop=ymax, length=nye))
    ze = collect(range(zmin, stop=zmax, length=nze))
    dx = abs(xe[2] - xe[1])
    dy = abs(ye[2] - ye[1])
    dz = abs(ze[2] - ze[1])

    fk = [f(x, y, z) for x in xk, y in yk, z in zk]
    fe = [f(x, y, z) for x in xe, y in ye, z in ze]

    @time s  = interpolate(hermite_x, hermite_y, hermite_z, fk, xe, ye, ze)

    max_error = maximum(abs.(s .- fe))

    l2_error  = sqrt(trapz((s - fe).^2, dx, dy, dz))

    println(" Hermite 3d max error : $max_error ")
    println(" Hermite 3d L2  error : $l2_error  ")

    @test true

end
