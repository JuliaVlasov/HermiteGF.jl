"""
    Radial( nodes::NodesType, epsilon, gamma )

Radial inteprolation

"""
mutable struct Radial <: InterpolationType

    nx         :: Int64
    epsilon    :: Float64
    xk         :: Array{Float64,1}
    colloc_mat :: Array{Float64,2}

    function Radial( nodes::NodesType, epsilon )
        
        xk         = nodes.xk
        nx         = nodes.nx
        colloc_mat = evaluate_radial(xk, xk, nx, epsilon)
        new( nx, epsilon, xk, colloc_mat )
        
    end

end

"""
    Radial( xe )

returns evaluation matrix
"""
function (interp::Radial)( xe::Array{Float64,1} )

    nx      = interp.nx
    epsilon = interp.epsilon
    xk      = interp.xk
    evaluate_radial(xk, xe, nx, epsilon)

end
