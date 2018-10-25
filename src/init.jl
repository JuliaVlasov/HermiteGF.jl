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
# Include all necessary files for the interpolation. This has to be ran every time
# when Julia is started.
###
if !isdefined(:interpolate_1D)
  include("helper_functions.jl")
  include("ndgrid.jl")
  include("evaluate_tensor_product.jl")
  include("evaluate_tensor_product_parallel.jl")
  include("evaluate_hermite.jl")
  include("interpolation.jl")
else
  print("You already included all necessary functions! If you changed the code of some function, recompile this function separately!")
end
