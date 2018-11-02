"""
    Uniform( xmin, xmax, nx )

Uniform nodes

"""
struct Uniform <: NodesType

    nx    :: Int64
    xmin  :: Float64
    xmax  :: Float64
    xk    :: Vector{Float64}

    function Uniform( xmin, xmax, nx )

        xk  = zeros(Float64, nx)
        xk .= range(-1, stop=1, length=nx);
        xk .= ((xk .+ 1) * (xmax - xmin)) /2 .+ xmin
        new( nx, xmin, xmax, xk )

    end

end
