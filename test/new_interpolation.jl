abstract type NodesType end

struct Chebyshev <: NodesType

    nx    :: Int64
    xmin  :: Float64
    xmax  :: Float64
    xk    :: Vector{Float64}

    function Chebyshev( nx, xmin, xmax )

        xk = generate_chebyshev_nodes(nx, xmin, xmax)
        new( nx, xmin, xmax, xk )

    end

end

struct Uniform <: NodesType

    nx    :: Int64
    xmin  :: Float64
    xmax  :: Float64
    xk    :: Vector{Float64}

    function Chebyshev( nx, xmin, xmax )

        xk = generate_uniform_nodes(nx, xmin, xmax)
        new( nx, xmin, xmax, xk )

    end

end

abstract type InterpolationType end

mutable struct Hermite <: InterpolationType
	
    dims               :: Int64
    epsilon            :: Float64
    gamma              :: Float64
    nodes              :: NodesType

    collocation_matrix :: Array{Float64,2}
    evaluation_matrix  :: Array{Float64,2}

    function Hermite( nodes, epsilon, gamma )

	dims = 1
        xk   = nodes.xk
	nx   = nodes.nx
	collocation_matrix = evaluate_hermite(xk, nx, epsilon, gamma)
        evaluation_matrix  = evaluate_hermite(xe, nx, epsilon, gamma)
	new( dims, epsilon, gamma, collocation_matrix, evaluation_matrix )

    end

end

mutable struct Radial <: InterpolationType

    epsilon :: Float64
    gamma   :: Float64
    collocation_matrix :: Array{Float64,2}
    evaluation_matrix  :: Array{Float64,2}

    function Hermite( nx, epsilon, gamma )
	dims = 1
	collocation_matrix_x = (xk, nx, epsilon)     -> evaluate_radial(xk, nx, epsilon, gamma)
        evaluation_matrix_x  = (xk, xe, nx, epsilon) -> evaluate_radial(xe, nx, epsilon, gamma)
	new( dims, epsilon, gamma, collocation_matrix, evaluation_matrix )
    end

end


function interpolate!(nodes_type::NodesType, 
		      interpolation_type::InterpolationType,
		      ne::Int64 )

   ϵ        = interpolation_type.epsilon
   nx       = nodes_type.nx
   xk       = generate_nodes_x(nx)
   min_dist = minimum(abs.(xk[2:end] .- xk[1:end-1]))

   x_collocation = interpolation_type.collocation_matrix_x(xk, nx, ϵ)
   x_evaluation  = interpolation_type.evaluation_matrix_x(xk, xe, nx, ϵ)

   # Precompute 1D interpolation matrix
   x_all = x_evaluation / x_collocation

   # Compute the interpolated values
   evaluate_s_1D(X_all, f)

end

## Setup the interpolated function.
#h(x) = cos.(x.^2)
#
## Setup the bounds of the interpolation domain.
#xmin = -1
#xmax = 1
#
#f  = h(x_k)
#xe = range(xmin, stop=xmax, length=ne1)
#
#dx = abs(xe[2] - xe[1])
#l2_error = sqrt(trapz_1D((s - f_vals).^2, dx))
#l1_error = maximum(abs.(s .- f_vals))
#
#f_vals = h(xx_e)
