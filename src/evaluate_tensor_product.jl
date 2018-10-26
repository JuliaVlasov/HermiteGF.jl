"""
In this file the functions for the serial computing of the interpolated values 's'
in 1-4D are implemented. Here the tensor representation of 's' is used (see Sec. 4.1 of the paper).

Computing s = (W_all ⊗ Z_all ⊗ Y_all ⊗ X_all) vec(F) (see Sec. 4.1 of the paper)
    
"""
function evaluate_s_4D(X_all, Y_all, Z_all, W_all, F)

    # Extract the number of collocation and evaluation points in each dimension
    Nx = size(X_all, 2);
    Ne1 = size(X_all, 1);
    Ny = size(Y_all, 2);
    Ne2 = size(Y_all, 1);
    Nz = size(Z_all, 2);
    Ne3 = size(Z_all, 1);
    Nw = size(W_all, 2);
    Ne4 = size(W_all, 1);
    
    # Initialize the result tensor
    s = zeros(Ne1, Ne2, Ne3, Ne4);
    
    @inbounds for col_dim_4 = 1:Nw
      @inbounds for eval_dim_4 = 1:Ne4
        @inbounds for col_dim_3 = 1:Nz
          @inbounds for eval_dim_3 = 1:Ne3
            @inbounds for col_dim_2 = 1:Ny
              @inbounds for eval_dim_2 = 1:Ne2
                @inbounds for col_dim_1 = 1:Nx
                  @inbounds for eval_dim_1 = 1:Ne1
                    s[eval_dim_1, eval_dim_2, eval_dim_3, eval_dim_4] = ( 
    		s[eval_dim_1,eval_dim_2, eval_dim_3, eval_dim_4] 
    		+ X_all[eval_dim_1, col_dim_1] * Y_all[eval_dim_2, col_dim_2]
    		* Z_all[eval_dim_3, col_dim_3] * W_all[eval_dim_4, col_dim_4]
		* F[col_dim_1, col_dim_2, col_dim_3, col_dim_4])
                  end
                end
              end
            end
          end
        end
      end
    end
    s
end

"""

Extract the number of collocation and evaluation points in each dimension 
Computing `s = (Z_all ⊗ Y_all ⊗ X_all) vec(F)` (see Sec. 4.1 of the paper)

"""
function evaluate_s_3D(X_all, Y_all, Z_all, F)

   Nx  = size(X_all, 2);
   Ne1 = size(X_all, 1);
   Ny  = size(Y_all, 2);
   Ne2 = size(Y_all, 1);
   Nz  = size(Z_all, 2);
   Ne3 = size(Z_all, 1);
   
   # Initialize the result tensor
   s = zeros(Ne1, Ne2, Ne3);
   
   @inbounds for col_dim_3 = 1:Nz
     @inbounds for eval_dim_3 = 1:Ne3
       @inbounds for col_dim_2 = 1:Ny
         @inbounds for eval_dim_2 = 1:Ne2
           @inbounds for col_dim_1 = 1:Nx
             @inbounds for eval_dim_1 = 1:Ne1
               s[eval_dim_1, eval_dim_2, eval_dim_3] = ( 
   	    s[eval_dim_1,eval_dim_2, eval_dim_3] 
   	    + X_all[eval_dim_1, col_dim_1]
   	    * Y_all[eval_dim_2, col_dim_2]
   	    * Z_all[eval_dim_3 ,col_dim_3]
	    * F[col_dim_1, col_dim_2, col_dim_3] )
             end
           end
         end
       end
     end
   end
   s

end

"""
Extract the number of collocation and evaluation points in each dimension

Computing `s = (Y_all ⊗ X_all) vec(F)` (see Sec. 4.1 of the paper)
"""
function evaluate_s_2D(X_all, Y_all, F)

    Nx = size(X_all, 2);
    Ne1 = size(X_all, 1);
    Ny = size(Y_all, 2);
    Ne2 = size(Y_all, 1);

    # Initialize the result tensor
    s = zeros(Ne1, Ne2);

    @inbounds for col_dim_2 = 1:Ny
      @inbounds for eval_dim_2 = 1:Ne2
        @inbounds for col_dim_1 = 1:Nx
          @inbounds for eval_dim_1 = 1:Ne1
            s[eval_dim_1, eval_dim_2] = ( s[eval_dim_1,eval_dim_2] 
	    + X_all[eval_dim_1, col_dim_1]
	    * Y_all[eval_dim_2, col_dim_2]*F[col_dim_1, col_dim_2])
          end
        end
      end
    end
    s
end

"""
Extract the number of collocation and evaluation points

Computing s =  X_all * vec(F) (see Sec. 4.1 of the paper)
"""
function evaluate_s_1D(X_all, F)

  Nx  = size(X_all, 2);
  Ne1 = size(X_all, 1);

  # Initialize the result tensor
  s = zeros(Ne1);

  @inbounds for col_dim_1 = 1:Nx
    @inbounds for eval_dim_1 = 1:Ne1
      s[eval_dim_1] = s[eval_dim_1] + X_all[eval_dim_1, col_dim_1]*F[col_dim_1];
    end
  end

  s

end
