module HermiteGF

  export Hermite, Radial
  export Chebyshev, Uniform
  export interpolate

  " Interpolation type (Hermite or Radial) "
  abstract type InterpolationType end

  " Node positions (Uniform or Chebyshev) "
  abstract type NodesType end

  include("chebyshev.jl")
  include("uniform.jl")
  include("hermite.jl")
  include("radial.jl")
  include("interpolation.jl")

end
