using TensorOperations

"""
    evaluate_s(x, y, z, v, w, f)

Computing s = (w ⊗ v ⊗ z  ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)
    
"""
function evaluate_s(x, y, z, v, w, f)

    # Extract the number of collocation and evaluation points in each dimension
    ne1 = size(x)[1]
    ne2 = size(y)[1]
    ne3 = size(z)[1]
    ne4 = size(v)[1]
    ne5 = size(w)[1]
    
    # Initialize the result tensor
    s = zeros(ne1, ne2, ne3, ne4, ne5)
    
    @tensor begin
        s[e1,e2,e3,e4,e5] = s[e1,e2,e3,e4,e5] + x[e1,c1]*y[e2,c2]*z[e3,c3]*v[e4,c4]*w[e5,c5]*f[c1,c2,c3,c4,c5]
    end
    s
end

"""
    evaluate_s(x, y, z, w, f)

Computing s = (w ⊗ z ⊗ y ⊗ x) vec(f) (see Sec. 4.1 of the paper)
    
"""
function evaluate_s(x, y, z, w, f)

    # Extract the number of collocation and evaluation points in each dimension
    ne1 = size(x)[1]
    ne2 = size(y)[1]
    ne3 = size(z)[1]
    ne4 = size(w)[1]
    
    # Initialize the result tensor
    s = zeros(ne1, ne2, ne3, ne4)
    
    @tensor begin
        s[e1,e2,e3,e4] = s[e1,e2,e3,e4] + x[e1,c1]*y[e2,c2]*z[e3,c3]*w[e4,c4]*f[c1,c2,c3,c4]
    end
    s
end

"""
    evaluate_s(x, y, z, f)

Extract the number of collocation and evaluation points in each dimension 
Computing `s = (z ⊗ y ⊗ x) vec(f)` (see Sec. 4.1 of the paper)

"""
function evaluate_s(x, y, z, f)

    ne1 = size(x)[1]
    ne2 = size(y)[1]
    ne3 = size(z)[1]
   
   # Initialize the result tensor
   s = zeros(ne1,ne2,ne3)
   
   @tensor begin
       s[e1,e2,e3] = s[e1,e2, e3]+x[e1,c1]*y[e2,c2]*z[e3,c3]*f[c1,c2,c3]
   end
   s

end

"""

    evaluate_s(x, y, f)

Extract the number of collocation and evaluation points in each dimension

Computing `s = (y ⊗ x) vec(f)` (see Sec. 4.1 of the paper)
"""
function evaluate_s(x, y, f)

    ne1 = size(x)[1]
    ne2 = size(y)[1]

    # Initialize the result tensor
    s = zeros(ne1, ne2)

    @tensor begin 
        s[e1,e2] = s[e1,e2]+x[e1,c1]*y[e2,c2]*f[c1,c2]
    end
    s
end

"""
    evaluate_s(x, f)

Extract the number of collocation and evaluation points

Computing s =  x * vec(f) (see Sec. 4.1 of the paper)
"""
function evaluate_s(x, f)

  # Initialize the result tensor
  s = similar(x)

  @tensor begin 
      s[e1] = s[e1] + x[e1,c1]*f[c1]
  end

  s

end
