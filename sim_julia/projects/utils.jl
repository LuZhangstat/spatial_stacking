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

function ldiagmul!(Da::AbstractVector, B::AbstractMatrix)
    # B = diagonal(Da) *B
    m,n = size(B)
    @boundscheck m == length(Da) || throw(BoundsError())
    
    @inbounds for j in 1:n
        for i in 1:m
            B[i, j] *= Da[i]
        end
    end
    return B
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
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, nk_k_list,
        y, X, XTX_list, XTy_list, priors)
    
    ## compute expectation of response in fold k ##
    
    # preallocation and precomputation
    L_grid_deltasq = length(deltasq_grid);
    out_put = Array{Float64, 2}(undef, nk_k_list[k], L_grid_deltasq);
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
            cholesky!(Symmetric(chol_inv_M, :U));
            u[:] = y[CV_ind_ls[k]] /  deltasq_pick;
        else
            chol_inv_M[1:p, 1:p] = XTX_list[k] / deltasq_pick + priors["inv_V_β"];
            chol_inv_M[1:p, (p+1):end] = X[:, CV_ind_ls[k]] / deltasq_pick;
            chol_inv_M[(p+1):end, (p+1):end] = invR_nk; 
            plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
            #print("in the inner for loop \n");
            cholesky!(Symmetric(chol_inv_M, :U));
            #print("right after inner chol \n");
            u[:] = [priors["inv_V_μ_β"] + XTy_list[k] / deltasq_pick; 
                y[CV_ind_ls[k]] /  deltasq_pick];
        end

        ldiv!(UpperTriangular(chol_inv_M)', u);  # u2 = chol_inv_M.L \ u;
        
        ldiv!(UpperTriangular(chol_inv_M), u); # compute the posterior expectation of β and z 

        if p == 0
            out_put[:, i2] = R_k_nk * (invR_nk * u);
        else
            out_put[:, i2] = 
                (X[:, CV_ind_hold_ls[k]]' * u[1:p]) + R_k_nk * (invR_nk * u[(p + 1):end]);
        end
                 
    end
    return out_put
end

function stacking_prediction_LP(coords, nu_pick, phi_pick, deltasq_grid, 
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, nk_k_list,
        y, X, XTX_list, XTy_list, y_sq_sum_list, priors, J = 300)
    
    ## compute expected log predictive density of response in fold k ##

    # preallocation and precomputation
    L_grid_deltasq = length(deltasq_grid);
    out_put = Array{Float64, 2}(undef, nk_k_list[k], L_grid_deltasq);
    chol_inv_M = Array{Float64, 2}(undef, nk_list[k] + p, nk_list[k] + p);
    u = Array{Float64, 1}(undef, nk_list[k] + p);
    a_star = 0.0; b_star = 0.0;
    σ2_sam = Array{Float64, 1}(undef, J);
    γ_sam = Array{Float64, 2}(undef, J, nk_list[k] + p);
    M_r = Array{Float64, 2}(undef, J, nk_k_list[k]);
    M_r_invRR_store = Array{Float64, 2}(undef, nk_list[k], nk_k_list[k]);
    if p > 0
        M_r_Xβ_store = Array{Float64, 2}(undef, J, nk_k_list[k]);
    end
    
    invR_nk = compute_invR_nk(coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    R_k_nk = compute_R_k_nk(coords[:, CV_ind_hold_ls[k]], 
        coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    
    for i2 in 1:L_grid_deltasq
        deltasq_pick = deltasq_grid[i2];
        if p == 0
            chol_inv_M[:] = invR_nk; 
            plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
            cholesky!(Symmetric(chol_inv_M, :U));
            u[:] = y[CV_ind_ls[k]] /  deltasq_pick;
        else
            chol_inv_M[1:p, 1:p] = XTX_list[k] / deltasq_pick + priors["inv_V_β"];
            chol_inv_M[1:p, (p+1):end] = X[:, CV_ind_ls[k]] / deltasq_pick;
            chol_inv_M[(p+1):end, (p+1):end] = invR_nk; 
            plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
            cholesky!(Symmetric(chol_inv_M, :U));
            u[:] = [priors["inv_V_μ_β"] + XTy_list[k] / deltasq_pick; y[CV_ind_ls[k]] /  deltasq_pick];
        end

        ldiv!(UpperTriangular(chol_inv_M)', u);  # u2 = chol_inv_M.L \ u;

        ## Stacking based on log point-wise predictive density
        if p == 0
            b_star = priors["bσ"] + 0.5 * (y_sq_sum_list[k] / deltasq_pick - norm(u)^2);
        else
            b_star = priors["bσ"] + 0.5 * (y_sq_sum_list[k] / deltasq_pick + 
                dot(priors["μβ"], priors["inv_V_μ_β"]) - norm(u)^2);
        end
        a_star = priors["aσ"] + nk_list[k] / 2;

        ## generate posterior samples ##
        rand!(InverseGamma(a_star, b_star), σ2_sam);

        ## compute the expected response on unobserved locations ##
        rand!(Normal(), γ_sam)   # each row is a sample
        ldiagmul!(sqrt.(σ2_sam), γ_sam);
        add_to_row!(γ_sam, u);          
        rdiv!(γ_sam, UpperTriangular(chol_inv_M)');            # γ_sam' = backsolve(chol_inv_M, gamma.sam)

        # M_r the matrix of log ratios (lp), first, store the expected responses on the held out fold
        if p == 0
            mul!(M_r_invRR_store, invR_nk, R_k_nk');
            mul!(M_r, γ_sam, M_r_invRR_store);
        else
            mul!(M_r_invRR_store, invR_nk, R_k_nk');
            mul!(M_r, γ_sam[:, (p+1):end], M_r_invRR_store);
            mul!(M_r_Xβ_store, γ_sam[:, (1:p)], X[:, CV_ind_hold_ls[k]]);
            M_r += M_r_Xβ_store;
        end

        # the matrix of log ratios (lp)
        minus_to_row!(M_r, y[CV_ind_hold_ls[k]]);
        ldiagmul!(sqrt.(1 ./ (deltasq_pick .* σ2_sam)), M_r);
        square_and_plus!(M_r, log(2*pi));
        add_to_column!(M_r, log.(deltasq_pick .* σ2_sam));
        M_r .*= -0.5;
        out_put[:, i2] = log.(mean(exp.(M_r), dims = 1))[1,:];
    end          
    return out_put
end

## compute the stacking weights ##
function QP_stacking_weight(Y_hat::Matrix, y::Vector)
    m,n = size(Y_hat)
    @boundscheck m == length(y) || throw(BoundsError())
    
    w = Variable(n);
    problem = minimize(sumsquares(y - Y_hat * w)); # objective
    problem.constraints += sum(w) == 1; # constraint
    problem.constraints += w >= 0; # constraint
    solver = () -> Mosek.Optimizer(LOG=0)
    solve!(problem, solver);
    return w.value[:]
end

function stacking_weight(lpd_point::AbstractMatrix)
    
    lp_m = mean(lpd_point);
    lpd_point .-=lp_m; # rescale the log-density for numerical stability
    exp_lpd_point = exp.(lpd_point);
    w = Variable(size(lpd_point, 2));
    problem = maximize(sum(log(exp_lpd_point * w))); # objective
    problem.constraints += sum(w) == 1; # constraint
    problem.constraints += w >= 0; # constraint
    solver = () -> Mosek.Optimizer(LOG=0);
    solve!(problem, solver);
    return w.value[:]
end