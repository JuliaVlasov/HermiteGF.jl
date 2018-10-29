"""
    Hermite( nodes::NodesType, epsilon:Real, gamma:Real )

Hermite inteprolation

"""
mutable struct Hermite <: InterpolationType

    nx         :: Int64
    epsilon    :: Float64
    gamma      :: Float64
    xk         :: Vector{Array{Float64,1}}
    colloc_mat :: AbstractVecOrMat{Float64}
    
    function Hermite( nodes::Vector{NodesType}, epsilon::Real, gamma::Real )

        xk   = nodes.xk
        nx   = nodes.nx
        colloc_mat = evaluate_hermite(xk, nx, epsilon, gamma)
        new( nx, epsilon, gamma, xk, colloc_mat )

    end

end

"""
    Hermite( xe )

returns evaluation matrix
"""
function (interp::Hermite)( xe::Array{Float64,1} )

   nx      = interp.nx
   epsilon = interp.epsilon
   gamma   = interp.gamma

   evaluate_hermite( xe, nx, epsilon, gamma)

end
