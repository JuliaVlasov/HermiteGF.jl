module HermiteGF

  export Hermite, Radial
  export Chebyshev, Uniform
  export interpolate

  include("helper_functions.jl")
  include("ndgrid.jl")
  include("evaluate_tensor_product.jl")
  include("evaluate_hermite.jl")

  " Interpolation type (Hermite or Radial) "
  abstract type InterpolationType end

  " Node positions (Uniform or Chebyshev) "
  abstract type NodesType end

  include("chebyshev.jl")
  include("uniform.jl")
  include("hermite.jl")
  include("radial.jl")

  include("interpolation1d.jl")
  include("interpolation2d.jl")
  include("interpolation3d.jl")
  include("interpolation4d.jl")
  include("interpolation5d.jl")

end
