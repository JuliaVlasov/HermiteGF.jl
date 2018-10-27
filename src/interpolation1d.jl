"""

     interpolate_1D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)

Computes an interpolant of a 1D function via Hermite-tensor method.

Arguments:
 - function_name : the name of the tested function.  Options_available: "f_3"
 - nodes_types : an array of two strings with the spacing of nodes in corresponding directions.
               Options available: :Uniform, :Chebyshev
 - epsilon : a vector of the values of shape parameters in corresponding directions.
           The values cannot be larger than 1!
 - N : a vector with the numbers of basis functions in corresponding directions.
 - Ne : a vector with the numbers of evaluation points in corresponding directions.
 - interpolation_type - a string with the name of preferred basis.
                      Options available: :Hermite, :Radial
 gamma - a factor allowing to scale the domain and the function accordingly.

  © Anna Yurova, 2017

"""
function interpolate_1D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)

   # Set the nodes types in all directions.
   nodes_type_x = nodes_types[1]


   ## ----------------------- TEST CASES -------------------------------------
   if function_name == :f_3
     # Setup the interpolated function.
     f(x) = cos.(x.^2)

     # Setup the bounds of the interpolation domain.
     xmin = -1
     xmax = 1
   else
     @error ("Wrong name of the test case! Possible options: :f_3")
   end

   ## -------------- COLLOCATION AND EVALUATION MATRICES ---------------------
   # Setup the basis used for the interpolation.
   # Declare the functions for the generating of the collocation and evaluation
   # matrices according to the chosen basis.
   if interpolation_type == :Hermite
     generate_collocation_matrix_x = (xk, N, epsilon) -> evaluate_hermite(xk, N, epsilon, gamma)
     generate_evaluation_matrix_x = (xk, xe, N, epsilon) -> evaluate_hermite(xe, N, epsilon, gamma)
   elseif interpolation_type == :Radial
     generate_collocation_matrix_x = (xk, N, epsilon) -> evaluate_radial(xk, xk, N, epsilon)
     generate_evaluation_matrix_x = (xk, xe, N, epsilon) -> evaluate_radial(xk, xe, N, epsilon)
   else
     @error ("Wrong basis type! Possible options are: Radial, Hermite.")
   end


   ## -------------------- NODES TYPES ---------------------------------------
   # Declare the nodes generating functions in all directions according to the specified by the user.
   # X - direction
   if nodes_type_x == :Chebyshev
     generate_nodes_x = (Nx) -> generate_chebyshev_nodes(Nx, xmin, xmax)
   elseif nodes_type_x == :Uniform
     generate_nodes_x = (Nx) -> generate_uniform_nodes(Nx, xmin, xmax)
   else
     @error ("Wrong node type! Possible options: Chebyshev, Uniform")
   end

   #  ---------------- INTERPOLATION PARAMETERS ------------------------
   # Setup the value of epsilon
   ep1 = epsilon[1]

   # Setup the number of collocation points
   Nx = N[1]

   #Setup the number of evaluation points
   Ne1 = Ne[1]

   # Generate node points
   xk = generate_nodes_x(Nx)
   min_dist = minimum(abs.(xk[2:end] .- xk[1:end-1]))
   xx_k = xk

   # Generate evaluation points
   xe = range(xmin, stop=xmax, length=Ne1)
   xx_e = xe

   ## ------------------- ACTUAL INTERPOLATION -------------------------------

   # Calculate collocation and evaluation matrices
   X_collocation = generate_collocation_matrix_x(xk, Nx, ep1)
   X_evaluation = generate_evaluation_matrix_x(xk, xe, Nx, ep1)

   # Precompute 1D interpolation matrix
   X_all = X_evaluation/X_collocation

   # Evaluate the function f at collocation points and evaluation points
   F = f(xx_k)
   f_vals = f(xx_e)

   # Compute the interpolated values
   s = evaluate_s_1D(X_all, F)

   ## ----------------- OBSERVABLES COMPUTATION ------------------------------

   # Maximum error
   Hermite_max_error = maximum(abs.(s .- f_vals))

   # L2 error. Using the fact that the evaluation grid is uniform
   dx = abs(xe[2] - xe[1])
   Hermite_L2_error = sqrt(trapz_1D((s - f_vals).^2, dx))

   [Hermite_max_error, Hermite_L2_error]

end

