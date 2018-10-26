###
# In this file the functions for the parallel computing of the interpolated values 's'
# in 4 and 5D are implemented. Here the tensor representation of 's' is used  (see Sec. 4.1 of the paper)
###

"""
This function is what we call from the rest of the code.
It distributes the work among processes and collects the results.
"""
function evaluate_s_4D_exec!(X_all, Y_all, Z_all, W_all, F, s)
    @sync begin
        for p in procs(s)
            @show p
            @async remotecall_wait(evaluate_s_4D_parallel_chunk!, p, X_all, Y_all, Z_all, W_all, F, s)
        end
    end
    return s;
end

"""
This function is what we call from the rest of the code.
It distributes the work among processes and collects the results.
"""
function evaluate_s_5D_exec!(X_all, Y_all, Z_all, W_all, V_all, F, s)
    @sync begin
        for p in procs(s)
            @show p
            @async remotecall_wait(evaluate_s_5D_parallel_chunk!, p, X_all, Y_all, Z_all, W_all, V_all, F, s)
        end
    end
    return s;
end

"""
This function retuns the (irange,jrange) indexes assigned to a worker.
In the case of the HermiteGF code, q corresponds to the tensor 's' containing
the values of the interpolant at the evaluation points
"""
@everywhere function myrange_4D(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0
    end

    # Split the number of points in 4-th dimension into <number of threads>
    # equal chuncks
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,4),nchunks+1)]
    return splits[idx]+1:splits[idx+1]
end

"""
This function retuns the (irange,jrange) indexes assigned to a worker
In the case of the HermiteGF code, q corresponds to the tensor 's' containing
the values of the interpolant at the evaluation points
"""
@everywhere function myrange_5D(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0
    end

    # Split the number of points in 5-th dimension into <number of threads>
    # equal chuncks
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,5),nchunks+1)]
    return splits[idx]+1:splits[idx+1]
end


"""
This is the kernel where the actual multiplication happens.
Every thread takes care of the slices within 'irange' of the evaluation points in 4-th dimension
"""
@everywhere function evaluate_s_4D_parallel!(X_all, Y_all, Z_all, W_all, F, s, irange)

    #@show (irange)  # display what irange each worker is taking care of

    # Extract the number of collocation and evaluation points in each dimension
    Nx = size(X_all, 2);
    Ne1 = size(X_all, 1);
    Ny = size(Y_all, 2);
    Ne2 = size(Y_all, 1);
    Nz = size(Z_all, 2);
    Ne3 = size(Z_all, 1);
    Nw = size(W_all, 2);
    Ne4 = size(W_all, 1);

    # Computing s(:, :, :, irange) = (W_all(irange, :) \otimes Z_all \otimes Y_all \otimes X_all) vec(F) (see Sec. 4.1 of the paper)
    for eval_dim_4 in irange
      @inbounds for col_dim_4 = 1:Nw
        @inbounds for col_dim_3 = 1:Nz
          @inbounds for eval_dim_3 = 1:Ne3
            @inbounds for col_dim_2 = 1:Ny
              @inbounds for eval_dim_2 = 1:Ne2
                @inbounds for col_dim_1 = 1:Nx
                  @inbounds for eval_dim_1 = 1:Ne1
                    s[eval_dim_1, eval_dim_2, eval_dim_3, eval_dim_4] = 
		    s[eval_dim_1,eval_dim_2, eval_dim_3, eval_dim_4] 
		    + X_all[eval_dim_1, col_dim_1]
		    * Y_all[eval_dim_2, col_dim_2]
		    * Z_all[eval_dim_3, col_dim_3]
		    * W_all[eval_dim_4, col_dim_4]
		    * F[col_dim_1, col_dim_2, col_dim_3, col_dim_4];
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
This is the kernel where the actual multiplication happens.
Every thread takes care of the slices within the 'irange' of the evaluation points in 5-th dimension
"""
@everywhere function evaluate_s_5D_parallel!(X_all, Y_all, Z_all, W_all, V_all, F, s, irange)
  
    #@show (irange)  # display what irange each worker is taking care of

    # Extract the number of collocation and evaluation points in each dimension
    Nx = size(X_all, 2);
    Ne1 = size(X_all, 1);
    Ny = size(Y_all, 2);
    Ne2 = size(Y_all, 1);
    Nz = size(Z_all, 2);
    Ne3 = size(Z_all, 1);
    Nw = size(W_all, 2);
    Ne4 = size(W_all, 1);
    Nv = size(V_all, 2);
    Ne5 = size(V_all, 1);

    # Computing s(:, :, :, :, irange) = (V_all(irange, :) \otimes W_all \otimes Z_all \otimes Y_all \otimes X_all) vec(F) (see Sec. 4.1 of the paper)
    for eval_dim_5 in irange
     @inbounds for col_dim_5 = 1:Nv
      @inbounds for col_dim_4 = 1:Nw
       @inbounds for eval_dim_4 = 1:Ne4
        @inbounds for col_dim_3 = 1:Nz
          @inbounds for eval_dim_3 = 1:Ne3
            @inbounds for col_dim_2 = 1:Ny
              @inbounds for eval_dim_2 = 1:Ne2
                @inbounds for col_dim_1 = 1:Nx
                  @inbounds for eval_dim_1 = 1:Ne1
                    s[eval_dim_1, eval_dim_2, eval_dim_3, eval_dim_4, eval_dim_5] = 
		    s[eval_dim_1,eval_dim_2, eval_dim_3, eval_dim_4, eval_dim_5] 
		    + X_all[eval_dim_1, col_dim_1]
		    * Y_all[eval_dim_2, col_dim_2]
		    * Z_all[eval_dim_3, col_dim_3]
		    * W_all[eval_dim_4, col_dim_4]
		    * V_all[eval_dim_5, col_dim_5]
		    * F[col_dim_1, col_dim_2, col_dim_3, col_dim_4, col_dim_5];
                  end
                end
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

# This is an interface simplifying the call of the function from the rest of the code
@everywhere evaluate_s_4D_parallel_chunk!(X_all, Y_all, Z_all, W_all, F, s) = evaluate_s_4D_parallel!(X_all, Y_all, Z_all, W_all, F, s, myrange_4D(s))
@everywhere evaluate_s_5D_parallel_chunk!(X_all, Y_all, Z_all, W_all, V_all, F, s) = evaluate_s_5D_parallel!(X_all, Y_all, Z_all, W_all, V_all, F, s, myrange_5D(s))
