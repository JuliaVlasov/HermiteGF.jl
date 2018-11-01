var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#HermiteGF.jl-Documentation-1",
    "page": "Documentation",
    "title": "HermiteGF.jl Documentation",
    "category": "section",
    "text": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials"
},

{
    "location": "#Types-1",
    "page": "Documentation",
    "title": "Types",
    "category": "section",
    "text": "```@autodocs Modules = [HermiteGF] Order   = [:type]"
},

{
    "location": "functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#HermiteGF.interpolate-Tuple{HermiteGF.InterpolationType,Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "HermiteGF.interpolate",
    "category": "method",
    "text": " interpolate( interp :: InterpolationType, \n              f  :: Array{Float64,1},\n              xe :: Array{Float64,1})\n\nComputes interpolation of a 1D function via Hermite-tensor method.\n\nArguments:\n\ninterp : inteprolation type (Hermite or Radial)\nf      : 1d array containing values at interpolation nodes (Chebyshev or Uniform)\nxe     : vector of evaluation points \n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "functions/#HermiteGF.interpolate-Tuple{HermiteGF.InterpolationType,HermiteGF.InterpolationType,Array{Float64,2},Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "HermiteGF.interpolate",
    "category": "method",
    "text": " interpolate( interp_x :: InterpolationType, \n              interp_y :: InterpolationType, \n              f  :: Array{Float64,2},\n              xe :: Array{Float64,1},\n              ye :: Array{Float64,1})\n\nComputes interpolation of a 2D function via Hermite-tensor method.\n\nArguments:\n\ninterp : inteprolation type (Hermite or Radial)\nf      : 2d array containing values at interpolation nodes (Chebyshev or Uniform)\nxe,ye  : vector of evaluation points \n\n\n\n\n\n"
},

{
    "location": "functions/#HermiteGF.interpolate-Tuple{HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,Array{Float64,3},Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "HermiteGF.interpolate",
    "category": "method",
    "text": " interpolate( interp_x :: InterpolationType, \n              interp_y :: InterpolationType, \n              interp_z :: InterpolationType, \n              f  :: Array{Float64,3},\n              xe :: Array{Float64,1},\n              ye :: Array{Float64,1},\n              ze :: Array{Float64,1})\n\nComputes interpolation of a 3D function via Hermite-tensor method.\n\nArguments:\n\ninterp_x : inteprolation type (Hermite or Radial)\ninterp_y : inteprolation type (Hermite or Radial)\ninterp_z : inteprolation type (Hermite or Radial)\nf        : 3d array containing values at interpolation nodes (Chebyshev or Uniform)\n\nComputing s = ( z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "functions/#HermiteGF.interpolate-Tuple{HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,Array{Float64,4},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "HermiteGF.interpolate",
    "category": "method",
    "text": " interpolate( interp_x :: InterpolationType, \n              interp_y :: InterpolationType, \n              interp_z :: InterpolationType, \n              interp_v :: InterpolationType, \n              f  :: Array{Float64,4},\n              xe :: Array{Float64,1},\n              ye :: Array{Float64,1},\n              ze :: Array{Float64,1},\n              ve :: Array{Float64,1})\n\nComputes interpolation of a 4D function via Hermite-tensor method.\n\nArguments:\n\ninterp_x : inteprolation type (Hermite or Radial)\ninterp_y : inteprolation type (Hermite or Radial)\ninterp_z : inteprolation type (Hermite or Radial)\ninterp_v : inteprolation type (Hermite or Radial)\nf        : 4d array containing values at interpolation nodes (Chebyshev or Uniform)\n\nComputing s = (v ⊗ z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "functions/#HermiteGF.interpolate-Tuple{HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,HermiteGF.InterpolationType,Array{Float64,5},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "HermiteGF.interpolate",
    "category": "method",
    "text": " interpolate( interp_x :: InterpolationType, \n              interp_y :: InterpolationType, \n              interp_z :: InterpolationType, \n              interp_v :: InterpolationType, \n              interp_w :: InterpolationType, \n              f  :: Array{Float64,4},\n              xe :: Array{Float64,1},\n              ye :: Array{Float64,1},\n              ze :: Array{Float64,1},\n              ve :: Array{Float64,1},\n              we :: Array{Float64,1})\n\nComputes interpolation of a 4D function via Hermite-tensor method.\n\nArguments:\n\ninterp_x : interpolation type (Hermite or Radial)\ninterp_y : interpolation type (Hermite or Radial)\ninterp_z : interpolation type (Hermite or Radial)\ninterp_v : interpolation type (Hermite or Radial)\ninterp_w : interpolation type (Hermite or Radial)\nf        : 5d arrray containing values at interpolation nodes (Chebyshev or Uniform)\n\nExtract the number of collocation and evaluation points in each dimension  Computing s = (w ⊗ v ⊗ z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "functions/#HermiteGF.evaluate_hermite-NTuple{4,Any}",
    "page": "Functions",
    "title": "HermiteGF.evaluate_hermite",
    "category": "method",
    "text": "evaluate_hermite(xk, n, epsilon, gamma)\n\nThis function returns the matrix He of the values of the HermiteGF basis functions. The computation is done via three term recurrence for Hermite functions with an argument gamma*x and then appropriate exponential scaling.\n\nMore details can be found in Section 5.1 of the paper\n\nSTABLE EVALUATION OF GAUSSIAN RADIAL BASIS  FUNCTIONS USING HERMITE POLYNOMIALS  by Anna Yurova and Katharina Kormann.\n\nThis function should be used both for computation of the collocation and evaluation matrices.\n\n\n\n\n\n"
},

{
    "location": "functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [HermiteGF]\nOrder   = [:function]"
},

]}
