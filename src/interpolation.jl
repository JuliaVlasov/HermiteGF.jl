using TensorOperations

"""

     interpolate( interp :: InterpolationType, 
                  f  :: Array{Float64,1},
                  xe :: Array{Float64,1})

Computes interpolation of a 1D function via Hermite-tensor method.

Arguments:
 - interp : inteprolation type (Hermite or Radial)
 - f      : 1d array containing values at interpolation nodes (Chebyshev or Uniform)
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

"""
     interpolate( interp_x :: InterpolationType, 
                  interp_y :: InterpolationType, 
                  f  :: Array{Float64,2},
                  xe :: Array{Float64,1},
                  ye :: Array{Float64,1})

Computes interpolation of a 2D function via Hermite-tensor method.

Arguments:
 - interp : inteprolation type (Hermite or Radial)
 - f      : 2d array containing values at interpolation nodes (Chebyshev or Uniform)
 - xe,ye  : vector of evaluation points 
"""
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



"""
     interpolate( interp_x :: InterpolationType, 
                  interp_y :: InterpolationType, 
                  interp_z :: InterpolationType, 
                  f  :: Array{Float64,3},
                  xe :: Array{Float64,1},
                  ye :: Array{Float64,1},
                  ze :: Array{Float64,1})

Computes interpolation of a 3D function via Hermite-tensor method.

Arguments:
 - interp_x : inteprolation type (Hermite or Radial)
 - interp_y : inteprolation type (Hermite or Radial)
 - interp_z : inteprolation type (Hermite or Radial)
 - f        : 3d array containing values at interpolation nodes (Chebyshev or Uniform)

Computing s = ( z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)
    
"""
function interpolate(interp_x :: InterpolationType, 
                     interp_y :: InterpolationType, 
                     interp_z :: InterpolationType, 
                     f        :: Array{Float64,3}, 
		     xe       :: Array{Float64,1},
		     ye       :: Array{Float64,1},
		     ze       :: Array{Float64,1})

    x = interp_x(xe) 
    y = interp_y(ye) 
    z = interp_z(ze) 

    @tensor begin
        s[e1,e2,e3] := x[e1,c1]*y[e2,c2]*z[e3,c3]*f[c1,c2,c3]
    end
    s
end

"""
     interpolate( interp_x :: InterpolationType, 
                  interp_y :: InterpolationType, 
                  interp_z :: InterpolationType, 
                  interp_v :: InterpolationType, 
                  f  :: Array{Float64,4},
                  xe :: Array{Float64,1},
                  ye :: Array{Float64,1},
                  ze :: Array{Float64,1},
                  ve :: Array{Float64,1})

Computes interpolation of a 4D function via Hermite-tensor method.

Arguments:
 - interp_x : inteprolation type (Hermite or Radial)
 - interp_y : inteprolation type (Hermite or Radial)
 - interp_z : inteprolation type (Hermite or Radial)
 - interp_v : inteprolation type (Hermite or Radial)
 - f        : 4d array containing values at interpolation nodes (Chebyshev or Uniform)

Computing s = (v ⊗ z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)
    
"""
function interpolate(interp_x :: InterpolationType, 
                     interp_y :: InterpolationType, 
                     interp_z :: InterpolationType, 
                     interp_v :: InterpolationType, 
                     f        :: Array{Float64,4}, 
		     xe       :: Array{Float64,1},
		     ye       :: Array{Float64,1},
		     ze       :: Array{Float64,1},
		     ve       :: Array{Float64,1})

    x = interp_x(xe) 
    y = interp_y(ye) 
    z = interp_z(ze) 
    v = interp_v(ve) 

    @tensor begin
        s[e1,e2,e3,e4] := x[e1,c1]*y[e2,c2]*z[e3,c3]*v[e4,c4]*f[c1,c2,c3,c4]
    end
    s
end

"""
     interpolate( interp_x :: InterpolationType, 
                  interp_y :: InterpolationType, 
                  interp_z :: InterpolationType, 
                  interp_v :: InterpolationType, 
                  interp_w :: InterpolationType, 
                  f  :: Array{Float64,4},
                  xe :: Array{Float64,1},
                  ye :: Array{Float64,1},
                  ze :: Array{Float64,1},
                  ve :: Array{Float64,1},
                  we :: Array{Float64,1})

Computes interpolation of a 4D function via Hermite-tensor method.

Arguments:
 - interp_x : interpolation type (Hermite or Radial)
 - interp_y : interpolation type (Hermite or Radial)
 - interp_z : interpolation type (Hermite or Radial)
 - interp_v : interpolation type (Hermite or Radial)
 - interp_w : interpolation type (Hermite or Radial)
 - f        : 5d arrray containing values at interpolation nodes (Chebyshev or Uniform)

Extract the number of collocation and evaluation points in each dimension 
Computing s = (w ⊗ v ⊗ z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)

"""
function interpolate(interp_x :: InterpolationType, 
                     interp_y :: InterpolationType, 
                     interp_z :: InterpolationType, 
                     interp_v :: InterpolationType, 
                     interp_w :: InterpolationType, 
                     f        :: Array{Float64,5}, 
		     xe       :: Array{Float64,1},
		     ye       :: Array{Float64,1},
		     ze       :: Array{Float64,1},
		     ve       :: Array{Float64,1},
		     we       :: Array{Float64,1})

    x = interp_x(xe) 
    y = interp_y(ye) 
    z = interp_z(ze) 
    v = interp_v(ve) 
    w = interp_w(we) 

    @tensor begin
        s[e1,e2,e3,e4,e5] := x[e1,c1]*y[e2,c2]*z[e3,c3]*v[e4,c4]*w[e5,c5]*f[c1,c2,c3,c4,c5]
    end
    s

end
