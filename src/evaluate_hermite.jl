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

function evaluate_hermite(xk, N, epsilon, gamma)
  # This function returns the matrix He of the values of the HermiteGF basis
  # functions. The computation is done via three term recurrence for Hermite
  # functions with an argument gamma*x and then appropriate exponential scaling.
  # More details can be found in Section 5.1 of the paper
  # "STABLE EVALUATION OF GAUSSIAN RADIAL BASIS
  #  FUNCTIONS USING HERMITE POLYNOMIALS"
  #  by Anna Yurova and Katharina Kormann.

  # This function should be used both for computation of the
  # collocation and evaluation matrices.

  # Initialize the result matrix
  result = zeros(Int(length(xk)), N);

  # Write the values of the first two Hermite functions to initialize the
  # three-term recurrence
  result[:, 1] = exp(-0.5.*((gamma*xk).^2));
  result[:, 2] = gamma*xk.*sqrt(2).*exp(-0.5.*((gamma*xk).^2));

  # Three term recurrence for the Hermite functions with argument gamma*x
  for i = 3:N
    result[:, i] = sqrt(2.0/(i-1)).*(gamma*xk).*result[:, i-1] - sqrt((i-2.0)/(i-1)).*result[:, i-2];
  end

  # Scaling the Hermite functions with the exponential factor
  # exp(-epsilon^2 x^2 + (gamma x)^2/2)
  for i=1:N
      result[:, i] = (pi^(1/4))*result[:, i].*exp((xk.^2).*(gamma*gamma*0.5 - epsilon^2));
  end
  return result;
end
