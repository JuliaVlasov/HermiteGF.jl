"""

    evaluate_hermite(xk, N, epsilon, gamma)

This function returns the matrix He of the values of the HermiteGF basis
functions. The computation is done via three term recurrence for Hermite
functions with an argument gamma*x and then appropriate exponential scaling.

More details can be found in Section 5.1 of the paper

*STABLE EVALUATION OF GAUSSIAN RADIAL BASIS
 FUNCTIONS USING HERMITE POLYNOMIALS*
 by Anna Yurova and Katharina Kormann.

This function should be used both for computation of the
collocation and evaluation matrices.

"""
function evaluate_hermite(xk, N, epsilon, gamma)

  # Initialize the result matrix
  result = zeros(Int(length(xk)), N)

  # Write the values of the first two Hermite functions to initialize the
  # three-term recurrence
  result[:,1] = exp.(-0.5.*((gamma*xk).^2))
  result[:,2] = gamma*xk.*sqrt(2).*exp.(-0.5.*((gamma*xk).^2))

  # Three term recurrence for the Hermite functions with argument gamma*x
  for i = 3:N
    result[:,i] = sqrt(2.0/(i-1)) .* (gamma*xk) .* result[:,i-1] - sqrt((i-2)/(i-1)) .* result[:,i-2]
  end

  # Scaling the Hermite functions with the exponential factor
  # exp(-epsilon^2 x^2 + (gamma x)^2/2)
  
  for i=1:N
      result[:,i] = (pi^(1/4)) * result[:,i] .* exp.((xk.^2) .* (gamma*gamma*0.5 - epsilon^2))
  end

  result

end
