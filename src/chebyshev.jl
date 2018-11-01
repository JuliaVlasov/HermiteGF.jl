"""
    Chebyshev( xmin, xmax, nx )

Chebyshev nodes

"""
struct Chebyshev <: NodesType

    nx    :: Int64
    xmin  :: Float64
    xmax  :: Float64
    xk    :: Vector{Float64}

    function Chebyshev( xmin, xmax, nx )

        xk  = zeros(Float64, nx)
        θ = (π:-π/(nx-1):-1e-14)
        xk  .= ((cos.(θ).+1)*(xmax-xmin))/2 .+ xmin
        new( nx, xmin, xmax, xk )

    end

end

