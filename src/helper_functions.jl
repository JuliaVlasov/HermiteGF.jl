###
# In this file the some helper functions are implemented. In particular,
# - evaluation of the radial basis
# - generation of chebyshev and uniform nodes (1D)
# - trapezoidal rule 1-5D for the UNIFORM grid
###

" Evaluate radial basis "
function evaluate_radial(vc, ve, N, epsilon)

    result = zeros(Int(length(ve)), Int(length(vc)));

    for i = 1:Int(length(ve)), i_c = 1:Int(length(vc)) # loop over center points

        result[i, i_c] = exp((-epsilon^2).*(ve[i] - vc[i_c]).^2);

    end

    result
end

" Generate chebyshev nodes on the interval [xmin xmax] "
function generate_chebyshev_nodes(N, xmin, xmax)

    theta1 = (pi:-pi/(N-1):-1e-14);
    xk = ((cos.(theta1).+1)*(xmax-xmin))/2 .+ xmin
    xk

end

" Generate uniform nodes on the interval [xmin xmax] "
function generate_uniform_nodes(n, xmin, xmax)

    nodes = range(-1, stop=1, length=n);
    xk = ((nodes .+ 1) * (xmax - xmin)) /2 .+ xmin
    xk

end

" Trapezoidal rule for the uniform grid for 1-5D "
function trapz_1D(values, dx)

    (sum(values) - values[1]/2 - values[end]/2)*dx

end

function trapz_2D(values, dx, dy)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
      temp[i] = trapz_1D(values[i, :], dy)
    end
    trapz_1D(temp, dx)

end

function trapz_3D(values, dx, dy, dz)

    temp = zeros(size(values, 1));
    for i = 1:size(values, 1)
        temp[i] = trapz_2D(values[i, :, :], dy, dz)
    end
    trapz_1D(temp, dx);

end

function trapz_4D(values, dx, dy, dz, dw)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
        temp[i] = trapz_3D(values[i, :, :, :], dy, dz, dw)
    end
    trapz_1D(temp, dx)

end

function trapz_5D(values, dx, dy, dz, dw, dv)

    temp = zeros(size(values, 1))
    for i = 1:size(values, 1)
      temp[i] = trapz_4D(values[i, :, :, :, :], dy, dz, dw, dv)
    end
    trapz_1D(temp, dx)

end