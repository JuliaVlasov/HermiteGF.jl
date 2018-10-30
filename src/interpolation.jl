"""

     interpolate( interp :: InterpolationType, 
                  f  :: Array{Float64,1},
                  xe :: Array{Float64,1})

Computes interpolation of a 1D function via Hermite-tensor method.

Arguments:
 - interp : inteprolation type (Hermite or Radial)
 - f      : 2d array containing values at interpolation nodes (Chebyshev or Uniform)
 - xe     : vector of evaluation points 

  © Anna Yurova, 2017

"""
function interpolate(interp :: InterpolationType, 
                     f      :: Array{Float64,1}, 
		     xe     :: Array{Float64,1} )

   x = interp(xe) 

   @tensor begin
      s[e] := x[e, c]*f[c]
   end
    
   s

end

function interpolate(interp_x :: InterpolationType, 
                     interp_y :: InterpolationType, 
                     f        :: Array{Float64,2}, 
		     xe       :: Array{Float64,1},
		     ye       :: Array{Float64,1})

   x = interp_x(xe) 
   y = interp_y(ye) 

   @tensor begin 
       s[e1,e2] := x[e1,c1]*y[e2,c2]*f[c1,c2]
   end
    
   s

end
