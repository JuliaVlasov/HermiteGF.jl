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
# In this file the test for the error scaling with the dimension is implemented.
# In particular, the dependence of the L2 error from the number of points per dimension.
###
using PyPlot

Nvec = 5:35; # Number of collocation point per dimension
ep = 0.1; # Value of the shape parameter (in this test the same in all directions, but can be different)
Ne = 53; # Number of evaluation points per dimension
gamma = 3; # scaling parameter gamma

# Initializing error vectors
# first row - maximum error
# second row - L2 error
errors_1D = zeros(2, length(Nvec));
errors_2D = zeros(2, length(Nvec));
errors_3D = zeros(2, length(Nvec));
errors_4D = zeros(2, length(Nvec));
errors_5D = zeros(2, length(Nvec));

# Loop through the different values of N
for ind = 1:length(Nvec)
  i = Nvec[ind];
  print("Current size: ", i, "\n");
  #=
  errors_5D[:, ind] = interpolate_5D("f_3", ["Chebyshev","Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep ep], [i i i i i], [Ne Ne Ne Ne Ne], "Hermite", gamma);
  errors_4D[:, ind] = interpolate_4D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep], [i i i i], [Ne Ne Ne Ne], "Hermite", gamma);=#
  errors_3D[:, ind] = interpolate_3D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep], [i i i], [Ne Ne Ne], "Hermite", gamma);
  errors_2D[:, ind] = interpolate_2D("f_3", ["Chebyshev", "Chebyshev"], [ep ep], [i i], [Ne Ne], "Hermite", gamma);
  errors_1D[:, ind] = interpolate_1D("f_3", ["Chebyshev"], [ep], [i], [Ne], "Hermite", gamma);
end

#=fig1 = figure("Maximum error")
#plot(Nvec, errors_4D[1, :], marker="o", label="4D");
plot(Nvec, errors_3D[1, :], marker="<", label="3D");
plot(Nvec, errors_2D[1, :], marker="<", label="2D");
plot(Nvec, errors_1D[1, :], marker=">", label="1D");
ax = gca() # Get the handle of the current axis
ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
legend(loc="upper right",fancybox="true")
xlabel("N")
ylabel("Maximum error")=#

# Plot the L2 error
fig2 = figure("L2 error scaling");
#plot(Nvec, errors_4D[2, :], marker="o", label="4D");
plot(Nvec, errors_3D[2, :], marker="<", label="3D");
plot(Nvec, errors_2D[2, :], marker="<", label="2D");
plot(Nvec, errors_1D[2, :], marker=">", label="1D");
ax = gca() # Get the handle of the current axis
ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
legend(loc="upper right",fancybox="true")
xlabel("N")
ylabel("L2 error")
