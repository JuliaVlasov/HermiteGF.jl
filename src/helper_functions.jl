#=The BSD license:
——————————————
Copyright (c) 2017 Anna Yurova, Max-Planck-Institut für Plasmaphysik.
All rights reserved.
Redistribution and use in source and binary forms, with or without modifi-
cation, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
* Neither the name of the Max-Planck-Institut für Plasmaphysik nor the
names of its contributors may be used to endorse or promote products derived
from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ”AS
IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABIL-
ITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL MAX-PLANCK-INSTITUT FÜR PLASMAPHYSIK
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EX-
EMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLI-
GENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
——————————————
In addition, we ask you to cite the following reference in scientific publica-
tions which contain results obtained with this software and developments:
A. Yurova, K. Kormann
“Stable evaluation of Gaussian radial basis functions using Hermite polyno-
mials”
=#

###
# In this file the some helper functions are implemented. In particular,
# - evaluation of the radial basis
# - generation of chebyshev and uniform nodes (1D)
# - trapezoidal rule 1-5D for the UNIFORM grid
###

function evaluate_radial(vc, ve, N, epsilon)
# Evaluate radial basis
result = zeros(Int(length(ve)), Int(length(vc)));
for i = 1:Int(length(ve)) # loop over evaluation points
    for i_c = 1:Int(length(vc)) # loop over center points
        result[i, i_c] = exp((-epsilon^2).*(ve[i] - vc[i_c]).^2);
    end
end
return result;
end

function generate_chebyshev_nodes(N, xmin, xmax)
# Generate chebyshev nodes on the interval [xmin xmax]
theta1 = (pi:-pi/(N-1):-1e-14);
xk = ((cos(theta1)+1)*(xmax-xmin))/2 + xmin;
return xk
end

function generate_uniform_nodes(N, xmin, xmax)
# Generate uniform nodes on the interval [xmin xmax]
nodes = linspace(-1, 1, N);
xk = ((nodes + 1) * (xmax - xmin)) /2 + xmin;
return xk
end

# Trapezoidal rule for the uniform grid for 1-5D
function trapz_1D(values, dx)
  return (sum(values) - values[1]/2 - values[end]/2)*dx;
end

function trapz_2D(values, dx, dy)
  temp = zeros(size(values, 1));
  for i = 1:size(values, 1)
    temp[i] = trapz_1D(values[i, :], dy);
  end
  return trapz_1D(temp, dx);
end

function trapz_3D(values, dx, dy, dz)
  temp = zeros(size(values, 1));
  for i = 1:size(values, 1)
    temp[i] = trapz_2D(values[i, :, :], dy, dz);
  end
  return trapz_1D(temp, dx);
end

function trapz_4D(values, dx, dy, dz, dw)
  temp = zeros(size(values, 1));
  for i = 1:size(values, 1)
    temp[i] = trapz_3D(values[i, :, :, :], dy, dz, dw);
  end
  return trapz_1D(temp, dx);
end

function trapz_5D(values, dx, dy, dz, dw, dv)
  temp = zeros(size(values, 1));
  for i = 1:size(values, 1)
    temp[i] = trapz_4D(values[i, :, :, :, :], dy, dz, dw, dv);
  end
  return trapz_1D(temp, dx);
end
