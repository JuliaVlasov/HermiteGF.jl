
"""

    interpolate_4D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)

Computes an interpolant of a 4D function via Hermite-tensor method.

Arguments:
 - function_name : the name of the tested function.
   Options_available: "f_3"
 - nodes_types : an array of three strings with the spacing of nodes in corresponding directions.
   Options available: :Uniform, :Chebyshev
 - epsilon : a vector of the values of shape parameters in corresponding directions.
           The values cannot be larger than 1!
 - N : a vector with the numbers of basis functions in corresponding directions.
 - Ne : a vector with the numbers of evaluation points in corresponding directions.
 - interpolation_type : a string with the name of preferred basis.
   Options available: :Hermite, :Radial
 gamma - a factor allowing to scale the domain and the function accordingly.

   © Anna Yurova, 2017
"""
function interpolate_4D(function_name, nodes_types, epsilon, N, Ne, interpolation_type, gamma)

  # Set the nodes types in all directions.
  nodes_type_x = nodes_types[1]
  nodes_type_y = nodes_types[2]
  nodes_type_z = nodes_types[3]
  nodes_type_w = nodes_types[4]


  ## ----------------------- TEST CASES -------------------------------------
  if function_name == :f_3

    # Setup the interpolated function.
    f = (x, y, z, w) -> cos(x.^2 + y.^2 + z.^2 + w.^2)
    # Setup the bounds of the interpolation domain.
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -1
    zmax = 1
    wmin = -1
    wmax = 1
  else
    @error ("Wrong name of the test case! Possible options: f_3")
  end

  ## -------------- COLLOCATION AND EVALUATION MATRICES ---------------------
  # Setup the basis used for the interpolation.
  # Declare the functions for the generating of the collocation and evaluation
  # matrices according to the chosed basis.
  if interpolation_type == "Hermite"
    generate_collocation_matrix_x = (xk, N, epsilon) -> evaluate_hermite(xk, N, epsilon, gamma)
    generate_evaluation_matrix_x = (xk, xe, N, epsilon) -> evaluate_hermite(xe, N, epsilon, gamma)

    generate_collocation_matrix_y = (xk, N, epsilon) -> evaluate_hermite(xk, N, epsilon, gamma)
    generate_evaluation_matrix_y = (xk, xe, N, epsilon) -> evaluate_hermite(xe, N, epsilon, gamma)

    generate_collocation_matrix_z = (zk, N, epsilon) -> evaluate_hermite(zk, N, epsilon, gamma)
    generate_evaluation_matrix_z = (zk, ze, N, epsilon) -> evaluate_hermite(ze, N, epsilon, gamma)

    generate_collocation_matrix_w = (wk, N, epsilon) -> evaluate_hermite(wk, N, epsilon, gamma)
    generate_evaluation_matrix_w = (wk, we, N, epsilon) -> evaluate_hermite(we, N, epsilon, gamma)
  elseif interpolation_type == "Radial"
    generate_collocation_matrix_x = (xk, N, epsilon) -> evaluate_radial(xk, xk, N, epsilon)
    generate_evaluation_matrix_x = (xk, xe, N, epsilon) -> evaluate_radial(xk, xe, N, epsilon)

    generate_collocation_matrix_y = (xk, N, epsilon) -> evaluate_radial(xk, xk, N, epsilon)
    generate_evaluation_matrix_y = (xk, xe, N, epsilon) -> evaluate_radial(xk, xe, N, epsilon)

    generate_collocation_matrix_z = (zk, N, epsilon) -> evaluate_radial(zk, zk, N, epsilon)
    generate_evaluation_matrix_z = (zk, ze, N, epsilon) -> evaluate_radial(zk, ze, N, epsilon)

    generate_collocation_matrix_w = (wk, N, epsilon) -> evaluate_radial(wk, wk, N, epsilon)
    generate_evaluation_matrix_w = (wk, we, N, epsilon) -> evaluate_radial(wk, we, N, epsilon)
  else
    error("Wrong basis type! Possible options are: Radial, Hermite.")
  end

  ## -------------------- NODES TYPES ---------------------------------------
  # Declare the nodes generating functions in all directions according to the specified by the user.
  # X - direction
  if nodes_type_x == "Chebyshev"
    generate_nodes_x = (Nx) -> generate_chebyshev_nodes(Nx, xmin, xmax)
  elseif nodes_type_x == "Uniform"
    generate_nodes_x = (Nx) -> generate_uniform_nodes(Nx, xmin, xmax)
  else
    error("Wrong node type! Possible options: Chebyshev, Uniform")
  end

  # Y - direction
  if nodes_type_y == "Chebyshev"
    generate_nodes_y = (Ny) -> generate_chebyshev_nodes(Ny, ymin, ymax)
  elseif nodes_type_y == "Uniform"
    generate_nodes_y = (Ny) -> generate_uniform_nodes(Ny,  ymin, ymax)
  else
    error("Wrong node type! Possible options: Chebyshev, Uniform")
  end

  # Z - direction
  if nodes_type_z == "Chebyshev"
    generate_nodes_z = (Nz) -> generate_chebyshev_nodes(Nz, zmin, zmax)
  elseif nodes_type_z == "Uniform"
    generate_nodes_z = (Nz) -> generate_uniform_nodes(Nz,  zmin, zmax)
  else
    error("Wrong node type! Possible options: Chebyshev, Uniform")
  end

  # W - direction
  if nodes_type_w == "Chebyshev"
    generate_nodes_w = (Nw) -> generate_chebyshev_nodes(Nw, wmin, wmax)
  elseif nodes_type_w == "Uniform"
    generate_nodes_w = (Nw) -> generate_uniform_nodes(Nw,  wmin, wmax)
  else
    error("Wrong node type! Possible options: Chebyshev, Uniform")
  end

## ---------------- INTERPOLATION PARAMETERS ------------------------
  # Setup the value of epsilon
  ep1 = epsilon[1]
  ep2 = epsilon[2]
  ep3 = epsilon[3]
  ep4 = epsilon[4]

  # Setup the number of collocation points.
  Nx = N[1]
  Ny = N[2]
  Nz = N[3]
  Nw = N[4]

  #Setup the number of evaluation points.
  Ne1 = Ne[1]
  Ne2 = Ne[2]
  Ne3 = Ne[3]
  Ne4 = Ne[4]

  # Generate node points.
  xk = generate_nodes_x(Nx)
  yk = generate_nodes_y(Ny)
  zk = generate_nodes_z(Nz)
  wk = generate_nodes_w(Nw)
  grid_col = ndgrid(xk, yk, zk, wk)
  xx_k = grid_col[1]
  yy_k = grid_col[2]
  zz_k = grid_col[3]
  ww_k = grid_col[4]

  # Generate evaluation points
  xe = linspace(xmin, xmax, Ne1)
  ye = linspace(ymin, ymax, Ne2)
  ze = linspace(zmin, zmax, Ne3)
  we = linspace(wmin, wmax, Ne4)
  grid_eval = ndgrid(xe, ye, ze, we)
  xx_e = grid_eval[1]
  yy_e = grid_eval[2]
  zz_e = grid_eval[3]
  ww_e = grid_eval[4]

  ## ------------------- ACTUAL INTERPOLATION -------------------------------

  # Calculate collocation and evaluation matrices
  X_collocation = generate_collocation_matrix_x(xk, Nx, ep1)
  X_evaluation = generate_evaluation_matrix_x(xk, xe, Nx, ep1)

  Y_collocation = generate_collocation_matrix_y(yk, Ny, ep2)
  Y_evaluation = generate_evaluation_matrix_y(yk, ye, Ny, ep2)

  Z_collocation = generate_collocation_matrix_z(zk, Nz, ep3)
  Z_evaluation = generate_evaluation_matrix_y(zk, ze, Nz, ep3)

  W_collocation = generate_collocation_matrix_w(wk, Nw, ep4)
  W_evaluation = generate_evaluation_matrix_w(wk, we, Nw, ep4)

  # Precompute 1D interpolation matrices
  X_all = X_evaluation/X_collocation
  Y_all = Y_evaluation/Y_collocation
  Z_all = Z_evaluation/Z_collocation
  W_all = W_evaluation/W_collocation

  # Evaluate the function f at collocation points and evaluation points
  F = f(xx_k, yy_k, zz_k, ww_k)
  f_vals = f(xx_e, yy_e, zz_e, ww_e)

  # Compute the interpolated values (in parallel!)
  s = SharedArray(Float64, Ne1, Ne2, Ne3, Ne4)
  s = evaluate_s_4D_exec!(X_all, Y_all, Z_all, W_all, F, s)

  ## ----------------- OBSERVABLES COMPUTATION ------------------------------

  # Maximum error
  Hermite_max_error = maximum(abs.(s .- f_vals))

  # L2 error. Using the fact that the evaluation grid is uniform
  dx = abs(xe[2] - xe[1])
  dy = abs(ye[2] - ye[1])
  dz = abs(ze[2] - ze[1])
  dw = abs(we[2] - we[1])
  Hermite_L2_error = sqrt(trapz_4D((s - f_vals).^2, dx, dy, dz, dw))

  return [Hermite_max_error, Hermite_L2_error]
end

