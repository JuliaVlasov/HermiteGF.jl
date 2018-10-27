var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#HermiteGF.interpolate_1D-NTuple{7,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.interpolate_1D",
    "category": "method",
    "text": " interpolate_1D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)\n\nComputes an interpolant of a 1D function via Hermite-tensor method.\n\nArguments:\n\nfunctionname : the name of the tested function.  Optionsavailable: \"f_3\"\nnodes_types : an array of two strings with the spacing of nodes in corresponding directions.             Options available: :Uniform, :Chebyshev\nepsilon : a vector of the values of shape parameters in corresponding directions.         The values cannot be larger than 1!\nN : a vector with the numbers of basis functions in corresponding directions.\nNe : a vector with the numbers of evaluation points in corresponding directions.\ninterpolation_type - a string with the name of preferred basis.                    Options available: :Hermite, :Radial\n\ngamma - a factor allowing to scale the domain and the function accordingly.\n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.interpolate_2D-NTuple{7,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.interpolate_2D",
    "category": "method",
    "text": "interpolate_2D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)\n\nComputes an interpolant of a 2D function via Hermite-tensor method.\n\nArguments:\n\nfunctionname : the name of the tested function.               Optionsavailable: \"f_3\"\nnodes_types : an array of two strings with the spacing of nodes in corresponding directions.             Options available: \"Uniform\", \"Chebyshev\"\nepsilon : a vector of the values of shape parameters in corresponding directions.         The values cannot be larger than 1!\nN : a vector with the numbers of basis functions in corresponding directions.\nNe : a vector with the numbers of evaluation points in corresponding directions.\ninterpolation_type : a string with the name of preferred basis.                    Options available: \"Hermite\", \"Radial\"\n\ngamma - a factor allowing to scale the domain and the function accordingly.\n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_hermite-NTuple{4,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_hermite",
    "category": "method",
    "text": "This function returns the matrix He of the values of the HermiteGF basis functions. The computation is done via three term recurrence for Hermite functions with an argument gamma*x and then appropriate exponential scaling.\n\nMore details can be found in Section 5.1 of the paper\n\nSTABLE EVALUATION OF GAUSSIAN RADIAL BASIS  FUNCTIONS USING HERMITE POLYNOMIALS  by Anna Yurova and Katharina Kormann.\n\nThis function should be used both for computation of the collocation and evaluation matrices.\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_radial-NTuple{4,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_radial",
    "category": "method",
    "text": "Evaluate radial basis \n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_s_1D-Tuple{Any,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_s_1D",
    "category": "method",
    "text": "Extract the number of collocation and evaluation points\n\nComputing s =  X_all * vec(F) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_s_2D-Tuple{Any,Any,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_s_2D",
    "category": "method",
    "text": "Extract the number of collocation and evaluation points in each dimension\n\nComputing s = (Y_all ⊗ X_all) vec(F) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_s_3D-NTuple{4,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_s_3D",
    "category": "method",
    "text": "Extract the number of collocation and evaluation points in each dimension  Computing s = (Z_all ⊗ Y_all ⊗ X_all) vec(F) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.evaluate_s_4D-NTuple{5,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.evaluate_s_4D",
    "category": "method",
    "text": "In this file the functions for the serial computing of the interpolated values \'s\' in 1-4D are implemented. Here the tensor representation of \'s\' is used (see Sec. 4.1 of the paper).\n\nComputing s = (Wall ⊗ Zall ⊗ Yall ⊗ Xall) vec(F) (see Sec. 4.1 of the paper)\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.generate_chebyshev_nodes-Tuple{Any,Any,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.generate_chebyshev_nodes",
    "category": "method",
    "text": "Generate chebyshev nodes on the interval [xmin xmax] \n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.generate_uniform_nodes-Tuple{Any,Any,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.generate_uniform_nodes",
    "category": "method",
    "text": "Generate uniform nodes on the interval [xmin xmax] \n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.interpolate_3D-NTuple{7,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.interpolate_3D",
    "category": "method",
    "text": "interpolate_3D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)\n\nComputes an interpolant of a 3D function via Hermite-tensor method.\n\nArguments:\n\nfunctionname : the name of the tested function.             Optionsavailable: \"f_3\"\nnodes_types : an array of three strings with the spacing of nodes in corresponding directions.           Options available: \"Uniform\", \"Chebyshev\"\nepsilon : a vector of the values of shape parameters in corresponding directions.       The values cannot be larger than 1!\nN : a vector with the numbers of basis functions in corresponding directions.\nNe : a vector with the numbers of evaluation points in corresponding directions.\ninterpolation_type : a string with the name of preferred basis.                  Options available: :Hermite, :Radial\ngamma : a factor allowing to scale the domain and the function accordingly.\n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.interpolate_4D-NTuple{7,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.interpolate_4D",
    "category": "method",
    "text": "interpolate_4D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)\n\nComputes an interpolant of a 4D function via Hermite-tensor method.\n\nArguments:\n\nfunctionname : the name of the tested function. Optionsavailable: \"f_3\"\nnodes_types : an array of three strings with the spacing of nodes in corresponding directions. Options available: :Uniform, :Chebyshev\nepsilon : a vector of the values of shape parameters in corresponding directions.         The values cannot be larger than 1!\nN : a vector with the numbers of basis functions in corresponding directions.\nNe : a vector with the numbers of evaluation points in corresponding directions.\ninterpolation_type : a string with the name of preferred basis. Options available: :Hermite, :Radial\n\ngamma - a factor allowing to scale the domain and the function accordingly.\n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.interpolate_5D-NTuple{7,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.interpolate_5D",
    "category": "method",
    "text": "interpolate_5D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)\n\nComputes an interpolant of a 5D function via Hermite-tensor method.\n\nArguments:\n\nfunctionname : the name of the tested function.  Optionsavailable: \"f_3\"\nnodes_types : an array of three strings with the spacing of nodes in corresponding directions.            Options available: :Uniform, :Chebyshev\nepsilon : a vector of the values of shape parameters in corresponding directions.        The values cannot be larger than 1!\nN : a vector with the numbers of basis functions in corresponding directions.\nNe : a vector with the numbers of evaluation points in corresponding directions.\ninterpolation_type : a string with the name of preferred basis. Options available: :Hermite, :Radial\ngamma - a factor allowing to scale the domain and the function accordingly.\n\n© Anna Yurova, 2017\n\n\n\n\n\n"
},

{
    "location": "index.html#HermiteGF.trapz_1D-Tuple{Any,Any}",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "HermiteGF.trapz_1D",
    "category": "method",
    "text": "Trapezoidal rule for the uniform grid for 1-5D \n\n\n\n\n\n"
},

{
    "location": "index.html#Stable-evaluation-of-Gaussian-radial-basis-functions-using-Hermite-polynomials-1",
    "page": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "title": "Stable evaluation of Gaussian radial basis functions using Hermite polynomials",
    "category": "section",
    "text": "Modules = [HermiteGF]\nOrder   = [:function, :type]interpolation_1D\ninterpolation_2D\ninterpolation_3D\ninterpolation_4D\ninterpolation_5D"
},

]}
