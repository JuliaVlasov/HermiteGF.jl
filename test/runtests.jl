###
# In this file the test for the error scaling with the dimension is implemented.
# In particular, the dependence of the L2 error from the number of points per dimension.
###
using Plots
pyplot()

using Test
# Number of collocation point per dimension
Nvec = 5:35 
# Value of the shape parameter (in this test the same in all directions, but can be different)
ep = 0.1
# Number of evaluation points per dimension
Ne = 53 
# scaling parameter gamma
gamma = 3 

# Initializing error vectors
# first row - maximum error
# second row - L2 error
errors_1D = zeros(2, length(Nvec));
errors_2D = zeros(2, length(Nvec));
errors_3D = zeros(2, length(Nvec));
errors_4D = zeros(2, length(Nvec));
errors_5D = zeros(2, length(Nvec));

# Loop through the different values of N
for (ind,i) in enumerate(Nvec)
  println("Current size: $i ");
  #=
  errors_5D[:, ind] = interpolate_5D("f_3", ["Chebyshev","Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep ep], [i i i i i], [Ne Ne Ne Ne Ne], "Hermite", gamma);
  errors_4D[:, ind] = interpolate_4D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep ep], [i i i i], [Ne Ne Ne Ne], "Hermite", gamma);=#
  #errors_3D[:, ind] = interpolate_3D("f_3", ["Chebyshev", "Chebyshev", "Chebyshev"], [ep ep ep], [i i i], [Ne Ne Ne], "Hermite", gamma);
  #errors_2D[:, ind] = interpolate_2D("f_3", ["Chebyshev", "Chebyshev"], [ep ep], [i i], [Ne Ne], "Hermite", gamma);
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
#plot(Nvec, errors_3D[2, :], marker="<", label="3D");
#plot(Nvec, errors_2D[2, :], marker="<", label="2D");
plot(Nvec, errors_1D[2, :], marker=">", label="1D");
ax = gca() # Get the handle of the current axis
ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
legend(loc="upper right",fancybox="true")
xlabel("N")
ylabel("L2 error")
