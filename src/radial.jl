"""
    Radial( nodes::NodesType, epsilon, gamma )

Radial interpolation

"""
mutable struct Radial <: InterpolationType

    nnodes     :: Int64
    epsilon    :: Float64
    nodes      :: Array{Float64,1}
    colloc_mat :: Array{Float64,2}

    function Radial( nodes_t::NodesType, epsilon::Float64 )
        
        nodes      = nodes_t.xk
        nnodes     = nodes_t.nx

        colloc_mat = exp.((-epsilon^2) * (nodes .- transpose(nodes)).^2)

        new( nnodes, epsilon, nodes, colloc_mat )
        
    end

end

"""
    Radial( xe )

returns evaluated function

"""
function (interp::Radial)( xe::Array{Float64,1} )

    epsilon = interp.epsilon
    xk      = interp.nodes

    (exp.((-epsilon^2) * (xe .- transpose(xk)).^2)) / interp.colloc_mat

end

