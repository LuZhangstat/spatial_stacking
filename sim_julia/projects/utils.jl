using SpecialFunctions

function Maternlu(x; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    
    # Matern function
    
    if x == 0.0
        σ2
    else
        a = ϕ * x;
        σ2 * 2^(1 - ν) / gamma(ν) * (a)^ν * besselk(ν, a)
    end
end

function MaternluU!(A::AbstractMatrix; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    
    # Matern function
    m,n = size(A)
    @boundscheck m == n || throw(BoundsError())
    a = 0.0;
    
    @inbounds for j in 1:n
        for i in 1:j
            if A[i,j] == 0.0
                A[i,j] = σ2
            else
                a = A[i,j] * ϕ;
                A[i,j] = σ2 * 2^(1 - ν) / gamma(ν) * (a)^ν * besselk(ν, a)
            end
        end
    end
    return A
end

function MaternluM!(A::AbstractMatrix; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    
    # Matern function
    m,n = size(A)
    a = 0.0;
    
    @inbounds for j in 1:n
        for i in 1:m
            if A[i,j] == 0.0
                A[i,j] = σ2
            else
                a = A[i,j] * ϕ;
                A[i,j] = σ2 * 2^(1 - ν) / gamma(ν) * (a)^ν * besselk(ν, a)
            end
        end
    end
    return A
end


function add_to_column!(A::AbstractMatrix, b::AbstractVector)
           m,n = size(A)
           @boundscheck m == length(b) || throw(BoundsError())
           
           @inbounds for j in 1:n
               for i in 1:m
                   A[i,j] += b[i]
               end
           end
           return A
       end

function minus_to_column!(A::AbstractMatrix,b::AbstractVector)
           m,n = size(A)
           @boundscheck m == length(b) || throw(BoundsError())
           
           @inbounds for j in 1:n
               for i in 1:m
                   A[i,j] -= b[i]
               end
           end
           return A
       end


function add_to_row!(A::AbstractMatrix,b::AbstractVector)
           m,n = size(A)
           @boundscheck n == length(b) || throw(BoundsError())
           
           @inbounds for j in 1:n
               for i in 1:m
                   A[i,j] += b[j]
               end
           end
           return A
       end

function minus_to_row!(A::AbstractMatrix,b::AbstractVector)
           m,n = size(A)
           @boundscheck n == length(b) || throw(BoundsError())
           
           @inbounds for j in 1:n
               for i in 1:m
                   A[i,j] -= b[j]
               end
           end
           return A
       end


function square_and_plus!(A::AbstractMatrix,b::AbstractFloat)
           m,n = size(A)           
           @inbounds for j in 1:n
               for i in 1:m
                   A[i,j] = A[i,j]^2
                   A[i,j] += b
               end
           end
           return A
       end

function compute_invR_nk(coords; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    ## compute the inverse Matern covariance matrix ##
    M = pairwise(Euclidean(), coords, dims = 2);
    MaternluU!(M; ν = ν, ϕ = ϕ, σ2 = σ2);
    M = Symmetric(M, :U);
    LinearAlgebra.inv!(cholesky!(M)) 
end

function compute_R_k_nk(coords1, coords2; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    ## compute the Matern covariance matrix ##
    M = pairwise(Euclidean(), coords1, coords2, dims = 2);
    MaternluM!(M, ν = ν, ϕ = ϕ, σ2 = σ2)
end

function plus_cI!(A::AbstractMatrix, c, p)
    # A[(p+1):end, (p+1):end] += cI  
    m,n = size(A)
    @boundscheck m == n || throw(BoundsError())
    
    @inbounds for i in (1 + p):n
            A[i,i] += c
    end
    return A
end


function stacking_prediction_LSE(coords, nu_pick, phi_pick, deltasq_grid, 
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, 
        y, X, XTX, XTy, inv_V_β, inv_V_μ_β)
    
    ## compute expectation of response in fold k ##
    
    # preallocation and precomputation
    L_grid_deltasq = length(deltasq_grid);
    out_put = Array{Float64, 2}(undef, length(CV_ind_hold_ls[k]), L_grid_deltasq);
    chol_inv_M = Array{Float64, 2}(undef, nk_list[k] + p, nk_list[k] + p);
    u = Array{Float64, 1}(undef, nk_list[k] + p);
    
    invR_nk = compute_invR_nk(coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    R_k_nk = compute_R_k_nk(coords[:, CV_ind_hold_ls[k]], 
        coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    
    # for each candidate value deltasq #
    for i2 in 1:L_grid_deltasq
        deltasq_pick = deltasq_grid[i2];
        if p == 0
            chol_inv_M[:] = invR_nk; 
            plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
            cholesky!(Symmetric(chol_inv_M, :U))
            u[:] = y[CV_ind_ls[k]] /  deltasq_pick;
        else
            chol_inv_M[1:p, 1:p] = XTX_list[k] / deltasq_pick + inv_V_β;
            chol_inv_M[1:p, (p+1):end] = X[:, CV_ind_ls[k]] / deltasq_pick;
            chol_inv_M[(p+1):end, (p+1):end] = invR_nk; 
            plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
            cholesky!(Symmetric(chol_inv_M, :U))
            u[:] = [inv_V_μ_β + XTy_list[k] / deltasq_pick; y[CV_ind_ls[k]] /  deltasq_pick];
        end

        ldiv!(UpperTriangular(chol_inv_M)', u);  # u2 = chol_inv_M.L \ u;
        
        ldiv!(UpperTriangular(chol_inv_M), u); # compute the posterior expectation of β and z u3 = chol_inv_M.L' \ u2;

        if p == 0
            out_put[:, i2] = R_k_nk * (invR_nk * u);
        else
            out_put[:, i2] = 
                (X[:, CV_ind_hold_ls[k]]' * u[1:p]) + R_k_nk * (invR_nk * u[(p + 1):end]);
        end
                 
    end
    return out_put
end
