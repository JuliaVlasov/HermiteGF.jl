function interpolate(interp::InterpolationType, 
        f::Array{Float64,1}, xe::Array{Float64,1} )

   ϵ  = interp.epsilon
   nx = interp.nx
   xk = interp.xk

   # Precompute 1D interpolation matrix
   x_all = interp(xe) / interp.colloc_mat

   # Compute the interpolated values
   ne = Int(length(xe))

   # Initialize the result tensor
   s = similar(xe)

   @tensor begin
      s[e] = s[e] + x_all[e, c]*f[c]
   end
    
   s

end
