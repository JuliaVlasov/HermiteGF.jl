module HermiteGF

  using Distributed
  using SharedArrays

  export interpolate_1D
  export interpolate_2D

  include("helper_functions.jl")
  include("ndgrid.jl")
  include("evaluate_tensor_product.jl")
  #@everywhere include("evaluate_tensor_product_parallel.jl")
  include("evaluate_hermite.jl")
  include("interpolation1d.jl")
  include("interpolation2d.jl")
  include("interpolation3d.jl")
  include("interpolation4d.jl")
  include("interpolation5d.jl")

end
