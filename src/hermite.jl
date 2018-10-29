"""
    Hermite( nodes::NodesType, epsilon, gamma )

Hermite inteprolation

"""
mutable struct Hermite <: InterpolationType

    nx         :: Vector{Int64}
    epsilon    :: Vector{Float64}
    gamma      :: Vector{Float64}
    xk         :: Vector{Array{Float64,1}}
    colloc_mat :: AbstractVecOrMat{Float64}
    
    function Hermite( nodes::Vector{NodesType}, 
            epsilon::Vector{Float64}, 
            gamma::Vector{Float64} )

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
