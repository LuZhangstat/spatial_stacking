using SpecialFunctions
using Distributions
using JuMP, Ipopt

function pick_deltasq(E_sigmasq, E_tausq; b = 2, p_ls = [0.1, 0.25, 0.5, 0.75, 0.9])
    # Use expectation of sigmasq and tausq to select alpha and beta
    alpha = b / E_sigmasq + 1
    beta = b / E_tausq + 1
    
    # Calculate the quantiles of the beta distribution
    quantile_ls = quantile(Beta(alpha, beta), p_ls)
    deltasq_cand = quantile_ls ./ (1 .- quantile_ls)
    
    return deltasq_cand
end


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

function row_squaresum_and_plus!(a::AbstractVector, A::AbstractMatrix,b::AbstractFloat)
           m,n = size(A) 
           @inbounds for j in 1:n
               for i in 1:m
                   if j == 1
                       a[i] = A[i,j]^2;
                   else
                       a[i] += A[i,j]^2;
                   end
               end
           end
           for i in 1:m
               a[i] += b;
           end
           return a
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

### code refining ###
function compute_R(coords; ν=1.0, ϕ = 6.0, σ2 = 1.0)
    ## compute the inverse Matern covariance matrix ##
    M = pairwise(Euclidean(), coords, dims = 2);
    MaternluU!(M; ν = ν, ϕ = ϕ, σ2 = σ2);
    M = Symmetric(M, :U);
end

function compute_invR_nk_lit(R_full, CholR, k, nk_k_ind, nk_list)
    if k == 1
        Lnk =  cholesky(R_full[nk_k_ind[k+1]:(nk_k_ind[K_fold+1]-1), 
            nk_k_ind[k+1]:(nk_k_ind[K_fold+1]-1)]).U;
    elseif (k>1)&(k<K_fold)
        Lnk = Array{Float64, 2}(undef, nk_list[k], nk_list[k]);
        Lnk[1:(nk_k_ind[k]-1), 1:(nk_k_ind[k]-1)] .= CholR[1:(nk_k_ind[k]-1), 1:(nk_k_ind[k]-1)];
        Lnk[1:(nk_k_ind[k]-1), nk_k_ind[k]:nk_list[k]] .= 
            CholR[1:(nk_k_ind[k]-1), nk_k_ind[k+1]:(nk_k_ind[K_fold+1]-1)];
                
        # Extract the relevant submatrix of `CholR`
        submatrix = CholR[nk_k_ind[k]:(nk_k_ind[K_fold+1]-1), nk_k_ind[k+1]:(nk_k_ind[K_fold+1]-1)]

        # Compute the upper triangular Cholesky decomposition efficiently
        Lnk[nk_k_ind[k]:nk_list[k], nk_k_ind[k]:nk_list[k]] .= cholesky(submatrix' * submatrix).U;
        
        Lnk = UpperTriangular(Lnk);
    else
        Lnk = UpperTriangular(CholR[1:(nk_k_ind[K_fold]-1), 1:(nk_k_ind[K_fold]-1)]);
    end
    LinearAlgebra.inv!(Lnk);
    return Lnk*Lnk';
end

function stacking_prediction_LSE_all(R_full, CholR, deltasq_grid, L_grid_deltasq,
        K_fold, p, nk_list, nk_k_list, nk_k_ind, y_ord, X_ord, XTX_list, XTy_list, priors)
    
    ## compute expectation of response in fold k ##
    out_put = Array{Float64, 2}(undef, nk_k_ind[K_fold+1]-1, L_grid_deltasq);
    for k in 1:K_fold
    # preallocation and precomputation
    
        chol_inv_M = Array{Float64, 2}(undef, nk_list[k] + p, nk_list[k] + p);
        u = Array{Float64, 1}(undef, nk_list[k] + p);
    
        invR_nk = compute_invR_nk_lit(R_full, CholR, k, nk_k_ind, nk_list);
        nk_indices = setdiff(1:(nk_k_ind[K_fold+1]-1), nk_k_ind[k]:(nk_k_ind[k+1]-1));
        R_k_nk = R_full[nk_k_ind[k]:(nk_k_ind[k+1]-1), nk_indices];
    
        # for each candidate value deltasq #
        for i2 in 1:L_grid_deltasq
            deltasq_pick = deltasq_grid[i2];
            if p == 0
                chol_inv_M[:] = invR_nk; 
                plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
                cholesky!(Symmetric(chol_inv_M, :U));
                u[:] = y_ord[nk_indices] /  deltasq_pick;
            else
                chol_inv_M[1:p, 1:p] = XTX_list[k] / deltasq_pick + priors["inv_V_β"];
                    chol_inv_M[1:p, (p+1):end] = X_ord[:, nk_indices] / deltasq_pick;
                chol_inv_M[(p+1):end, (p+1):end] = invR_nk; 
                plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
                #print("in the inner for loop \n");
                cholesky!(Symmetric(chol_inv_M, :U));
                #print("right after inner chol \n");
                u[:] = [priors["inv_V_μ_β"] + XTy_list[k] / deltasq_pick; 
                    y_ord[nk_indices] /  deltasq_pick];
            end

            ldiv!(UpperTriangular(chol_inv_M)', u);  # u2 = chol_inv_M.L \ u;
        
            ldiv!(UpperTriangular(chol_inv_M), u); # compute the posterior expectation of β and z 

            if p == 0
                out_put[nk_k_ind[k]:(nk_k_ind[k+1]-1), i2] = R_k_nk * (invR_nk * u);
            else
                out_put[nk_k_ind[k]:(nk_k_ind[k+1]-1), i2] = 
                    (X_ord[:, nk_k_ind[k]:(nk_k_ind[k+1]-1)]' * u[1:p]) + R_k_nk * (invR_nk * u[(p + 1):end]);
            end     
        end
    end
    return out_put
end

function stacking_prediction_LP_all(R_full, CholR, deltasq_grid, L_grid_deltasq,
        K_fold, p, nk_list, nk_k_list, nk_k_ind, y_ord, X_ord, XTX_list, XTy_list, 
        y_sq_sum_list, priors)
    
    ## compute expected log predictive density of response in fold k ##

    # preallocation and precomputation
    out_put = Array{Float64, 2}(undef, nk_k_ind[K_fold+1]-1, L_grid_deltasq);
    
    for k in 1:K_fold
        chol_inv_M = Array{Float64, 2}(undef, nk_list[k] + p, nk_list[k] + p);
        u = Array{Float64, 1}(undef, nk_list[k] + p);
    
        invR_nk = compute_invR_nk_lit(R_full, CholR, k, nk_k_ind, nk_list);
        nk_indices = setdiff(1:(nk_k_ind[K_fold+1]-1), nk_k_ind[k]:(nk_k_ind[k+1]-1));
        R_k_nk = R_full[nk_k_ind[k]:(nk_k_ind[k+1]-1), nk_indices];
        
        a_star = 0.0; b_star = 0.0;
        
        
        H = Array{Float64, 2}(undef, nk_k_list[k], nk_list[k] + p);
        y_U_store = Array{Float64, 1}(undef, nk_k_list[k]);
        Vs_store = Array{Float64, 1}(undef, nk_k_list[k]);
        one_v = ones(nk_list[k] + p);
    
        for i2 in 1:L_grid_deltasq
            deltasq_pick = deltasq_grid[i2];
            if p == 0
                chol_inv_M[:] = invR_nk; 
                plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
                cholesky!(Symmetric(chol_inv_M, :U));
                u[:] = y_ord[nk_indices] /  deltasq_pick;
            else
                chol_inv_M[1:p, 1:p] = XTX_list[k] / deltasq_pick + priors["inv_V_β"];
                chol_inv_M[1:p, (p+1):end] = X_ord[:, nk_indices] / deltasq_pick;
                chol_inv_M[(p+1):end, (p+1):end] = invR_nk; 
                plus_cI!(chol_inv_M, 1 / deltasq_pick, p);
                cholesky!(Symmetric(chol_inv_M, :U));
                u[:] = [priors["inv_V_μ_β"] + XTy_list[k] / deltasq_pick; 
                    y_ord[nk_indices] /  deltasq_pick];
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
        
            ldiv!(UpperTriangular(chol_inv_M), u); # compute the posterior expectation of β and z 

            lp_c = -0.5 * log(2 * pi) + loggamma(a_star + 0.5) - 
                loggamma(a_star) + a_star * log(b_star); # compute the constant related to a_star b_star
            if p == 0
                mul!(H, R_k_nk, invR_nk);
            else
                H[:, 1:p] = X_ord[:, nk_k_ind[k]:(nk_k_ind[k+1]-1)]';
                H[:, (p+1):end] = R_k_nk * invR_nk;
            end

            # posterior expectations for observations in fold k 
            mul!(y_U_store, H, u); 
            # variance of the marignal posterior predictive distr
            rdiv!(H, UpperTriangular(chol_inv_M));
            row_squaresum_and_plus!(Vs_store, H, deltasq_pick);
            out_put[nk_k_ind[k]:(nk_k_ind[k+1]-1), i2] = lp_c .- 0.5 * log.(Vs_store) .- (a_star + 0.5) .* 
                log.(b_star .+ (y_ord[nk_k_ind[k]:(nk_k_ind[k+1]-1)] - y_U_store).^2 ./ (2 * Vs_store));
        end
    end
    return out_put
end

## compute the stacking weights ##
function QP_stacking_weight(Y_hat::Matrix, y::Vector)
    m, n = size(Y_hat);
    @boundscheck m == length(y) || throw(BoundsError());

    # Initialize the model
    model = Model(Mosek.Optimizer);
    JuMP.set_optimizer_attribute(model, "LOG", 0);  # Set the LOG attribute for Mosek

    # Define the variables with non-negativity constraints
    @variable(model, w[1:n] >= 0);

    # Define the objective
    @objective(model, Min, sum((y[i] - sum(Y_hat[i, j] * w[j] for j in 1:n))^2 for i in 1:m));

    # Add the constraint that weights sum to 1
    @constraint(model, sum(w) == 1);

    # Solve the problem
    optimize!(model);

    # Return the solution
    if termination_status(model) == MOI.OPTIMAL
        return value.(w)
    else
        error("Optimization problem did not solve successfully: $(termination_status(model))")
    end
end

function stacking_weight(lpd_point::AbstractMatrix)
    lp_m = mean(lpd_point);
    lpd_point .-= lp_m;  # Rescale the log-density for numerical stability
    exp_lpd_point = exp.(lpd_point);

    model = Model(Ipopt.Optimizer); #Mosek is deprecated. Use Ipopt (Interior Point OPTimizer), a software package for large-scale nonlinear optimization.
    JuMP.set_optimizer_attribute(model, "print_level", 0);  # Equivalent to Mosek's LOG=0 for Ipopt
    
    JuMP.register(model, :sum, 1, sum, autodiff = true);

    @variable(model, w[1:size(lpd_point, 2)] >= 0);  # Define variables with non-negativity constraint
    #@NLobjective(model, Max, sum(log.(exp_lpd_point * w)));  # Define the nonlinear objective
    @NLobjective(model, Max, sum(log(sum(exp_lpd_point[i, j] * w[j] 
        for j in 1:size(lpd_point, 2))) for i in 1:size(lpd_point, 1)));

    @constraint(model, sum(w) == 1);  # Add the constraints

    optimize!(model);  # Solve the optimization problem

    if termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL
        return value.(w);  # Return the solution
    else
        error("Optimization problem did not solve successfully: $(termination_status(model))");
    end
end

function sp_stacking_K_fold_V2(X, y, coords, deltasq_grid, phi_grid,
    nu_grid, Priors; K_fold = 10, seed = 123, label = "LSE")
    ## experimential multithread version of sp_stacking_K_fold
    
    Random.seed!(seed);
    # pre-computation and pre-allocation #
    N = size(X, 2);
    CV_ind_ls = collect(Kfold(N, K_fold)); # index of train data in CV
    CV_ind_hold_ls = [setdiff(1:N, CV_ind_ls[k]) for k in 1:K_fold]; # index of held-out data in CV
    N_grid = length(deltasq_grid) * length(phi_grid) * length(nu_grid);
    nk_list = [length(x) for x in CV_ind_ls]; # be careful, different from the nk_list in my R code
    nk_k_list = [(N - x) for x in nk_list];   # This is the nk_list in my R code
    nk_k_ind = vcat(1, [sum(nk_k_list[1:k])+1 for k in 1:K_fold]); # for partition the reorders data
    CV_reorder_ind = vcat([CV_ind_hold_ls[k] for k in 1:K_fold]...);

    if X == Nothing()
        p = 0;
    else
        p = size(X, 1);
        Priors["inv_V_μ_β"] = Priors["inv_V_β"] * Priors["μβ"];
        X_ord = X[:, CV_reorder_ind]; y_ord = y[CV_reorder_ind];
        XTX = X_ord * X_ord'; XTy = X_ord * y_ord;
        XTX_list = [XTX - X_ord[:, nk_k_ind[k]:(nk_k_ind[k+1]-1)] * X_ord[:, nk_k_ind[k]:(nk_k_ind[k+1]-1)]' 
            for k in 1:K_fold];
        XTy_list = [XTy - X_ord[:, nk_k_ind[k]:(nk_k_ind[k+1]-1)] * y_ord[nk_k_ind[k]:(nk_k_ind[k+1]-1)] 
            for k in 1:K_fold];
    end

    if label == "LSE"
        y_expect = Array{Float64, 2}(undef, N, N_grid);
    elseif label == "LP"
        lp_expect = Array{Float64, 2}(undef, N, N_grid);
        y_sq_sum_list = [norm(y_ord[setdiff(1:N, nk_k_ind[k]:(nk_k_ind[k+1]-1))])^2 for k in 1:K_fold];
    else 
        print("label has to be LSE or LP");
    end

    grid_phi_nu = vcat([[x y] for x in phi_grid, y in nu_grid]...);
    grid_all = vcat([[x y z] for x in phi_grid, y in nu_grid, z in deltasq_grid]...);
    L_grid_deltasq = length(deltasq_grid);
        
    ## Compute expectation for stacking ##
    #println("start: ");
    #prog = Progress(size(grid_phi_nu, 1), 1, "Computing initial pass...", 50)
    Threads.@threads for i1 in 1:size(grid_phi_nu, 1)
        #println("Thread $(Threads.threadid()) handles $i1");
        #t_start = time()
        phi_pick = grid_phi_nu[i1, 1];
        nu_pick = grid_phi_nu[i1, 2];
        coords_ord = coords[:, CV_reorder_ind]; # order by 10 folds
        R_full = compute_R(coords_ord, ν = nu_pick, ϕ = phi_pick);
        CholR = cholesky(R_full).U;
        if label == "LSE"
            y_expect[:, (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                    stacking_prediction_LSE_all(R_full, CholR, deltasq_grid, L_grid_deltasq,
                    K_fold, p, nk_list, nk_k_list, nk_k_ind, y_ord, X_ord, XTX_list, XTy_list, 
                    Priors); 
        else 
            lp_expect[:, (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                stacking_prediction_LP_all(R_full, CholR, deltasq_grid, L_grid_deltasq,
                    K_fold, p, nk_list, nk_k_list, nk_k_ind, y_ord, X_ord, XTX_list, XTy_list, 
                    y_sq_sum_list, Priors);
        end
        #println("Thread $(Threads.threadid()) finished $i1 in $(time() - t_start) seconds");
        #next!(prog)
    end
    
    # compute stacking weights
    if label == "LSE"
        w = QP_stacking_weight(y_expect, y_ord);
    else
        w = stacking_weight(lp_expect);
    end
    
    #return grid_all, weights, and label 
    return Dict(:grid_all => grid_all, :w => w, :label => label)
end

                        
                        

#---- Below is the old code that is not used in the time comparison ---#
#=                        
function stacking_prediction_LSE(coords, nu_pick, phi_pick, deltasq_grid, 
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, nk_k_list,
        y, X, XTX_list, XTy_list, priors)
    
    ## compute expectation of response in fold k ##
    
    # preallocation and precomputation
    #L_grid_deltasq = length(deltasq_grid);
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

function stacking_prediction_LP_MC(coords, nu_pick, phi_pick, deltasq_grid, 
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, nk_k_list,
        y, X, XTX_list, XTy_list, y_sq_sum_list, priors, J = 300)
    
    ## compute expected log predictive density of response in fold k, monte carlo version ##

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
    
    #Random.seed!(seed);
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


function stacking_prediction_LP(coords, nu_pick, phi_pick, deltasq_grid, 
        L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, nk_list, nk_k_list,
        y, X, XTX_list, XTy_list, y_sq_sum_list, priors, J = 300, MC = false)
    
    ## compute expected log predictive density of response in fold k ##

    # preallocation and precomputation
    L_grid_deltasq = length(deltasq_grid);
    out_put = Array{Float64, 2}(undef, nk_k_list[k], L_grid_deltasq);
    chol_inv_M = Array{Float64, 2}(undef, nk_list[k] + p, nk_list[k] + p);
    u = Array{Float64, 1}(undef, nk_list[k] + p);
    a_star = 0.0; b_star = 0.0;
    if MC == true
        σ2_sam = Array{Float64, 1}(undef, J);
        γ_sam = Array{Float64, 2}(undef, J, nk_list[k] + p);
        M_r = Array{Float64, 2}(undef, J, nk_k_list[k]);
        if p > 0
            M_r_Xβ_store = Array{Float64, 2}(undef, J, nk_k_list[k]);
        end
        M_r_invRR_store = Array{Float64, 2}(undef, nk_list[k], nk_k_list[k]);
    else
        H = Array{Float64, 2}(undef, nk_k_list[k], nk_list[k] + p);
        y_U_store = Array{Float64, 1}(undef, nk_k_list[k]);
        Vs_store = Array{Float64, 1}(undef, nk_k_list[k]);
        one_v = ones(nk_list[k] + p);
    end

    invR_nk = compute_invR_nk(coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    R_k_nk = compute_R_k_nk(coords[:, CV_ind_hold_ls[k]], 
        coords[:, CV_ind_ls[k]], ν = nu_pick, ϕ = phi_pick);
    
    #Random.seed!(seed);
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
        
        if MC == true
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
        else
            ldiv!(UpperTriangular(chol_inv_M), u); # compute the posterior expectation of β and z 

            lp_c = -0.5 * log(2 * pi) + loggamma(a_star + 0.5) - 
              loggamma(a_star) + a_star * log(b_star); # compute the constant related to a_star b_star
            if p == 0
                mul!(H, R_k_nk, invR_nk);
            else
                H[:, 1:p] = X[:, CV_ind_hold_ls[k]]';
                H[:, (p+1):end] = R_k_nk * invR_nk;
            end

            # posterior expectations for observations in fold k 
            mul!(y_U_store, H, u); 
            # variance of the marignal posterior predictive distr
            rdiv!(H, UpperTriangular(chol_inv_M));
            row_squaresum_and_plus!(Vs_store, H, deltasq_pick);
            out_put[:, i2] = lp_c .- 0.5 * log.(Vs_store) .- (a_star + 0.5) .* 
                log.(b_star .+ (y[CV_ind_hold_ls[k]] - y_U_store).^2 ./ (2 * Vs_store));
        end
    end          
    return out_put
end
=#


#=
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
=#


#=
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
=#


#=
function sp_stacking_K_fold(X, y, coords, deltasq_grid, phi_grid,
    nu_grid, Priors; K_fold = 10, seed = 123, label = "LSE", J = 300)
    
    Random.seed!(seed);
    # pre-computation and pre-allocation #
    N = size(X, 2);
    CV_ind_ls = collect(Kfold(N, K_fold)); # index of train data in CV
    CV_ind_hold_ls = [setdiff(1:N, CV_ind_ls[k]) for k in 1:K_fold]; # index of held-out data in CV
    N_grid = length(deltasq_grid) * length(phi_grid) * length(nu_grid);
    nk_list = [length(x) for x in CV_ind_ls]; # be careful, different from the nk_list in my R code
    nk_k_list = [(N - x) for x in nk_list];   # This is the nk_list in my R code


    if X == Nothing()
        p = 0;
    else
        p = size(X, 1);
        Priors["inv_V_μ_β"] = Priors["inv_V_β"] * Priors["μβ"];
        XTX = X * X'; XTy = X * y;
        XTX_list = [XTX - X[:, CV_ind_hold_ls[k]] * X[:, CV_ind_hold_ls[k]]' for k in 1:K_fold];
        XTy_list = [XTy - X[:, CV_ind_hold_ls[k]] * y[CV_ind_hold_ls[k]] for k in 1:K_fold];
    end

    if label == "LSE"
        y_expect = Array{Float64, 2}(undef, N, N_grid);
    elseif label == "LP"
        lp_expect = Array{Float64, 2}(undef, N, N_grid);
        y_sq_sum_list = [norm(y[CV_ind_ls[k]])^2 for k in 1:K_fold];
    else 
        print("label has to be LSE or LP");
    end

    grid_phi_nu = vcat([[x y] for x in phi_grid, y in nu_grid]...);
    grid_all = vcat([[x y z] for x in phi_grid, y in nu_grid, z in deltasq_grid]...);
    L_grid_deltasq  = length(deltasq_grid);
    
    ## Compute expectation for stacking ##
    prog = Progress(size(grid_phi_nu, 1), 1, "Computing initial pass...", 50)
    for i1 in 1:size(grid_phi_nu, 1)
        #println("Thread $(Threads.threadid()) handles $i1");
        #t_start = time()
        phi_pick = grid_phi_nu[i1, 1];
        nu_pick = grid_phi_nu[i1, 2];
        #Threads.@threads 
        for k in 1:K_fold
            Random.seed!(seed + (i1 - 1) * K_fold + k);
            if label == "LSE"
                y_expect[CV_ind_hold_ls[k], 
                    (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                stacking_prediction_LSE(coords, nu_pick, phi_pick, deltasq_grid, 
                    L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, 
                    nk_list, nk_k_list, y, X, XTX_list, XTy_list, Priors);
            else
                lp_expect[CV_ind_hold_ls[k], 
                    (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                stacking_prediction_LP(coords, nu_pick, phi_pick, deltasq_grid, 
                    L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, 
                    nk_list, nk_k_list, y, X, XTX_list, XTy_list,
                    y_sq_sum_list, Priors, J);
            end     
        end
        #println("Thread $(Threads.threadid()) finished $i1 in $(time() - t_start) seconds");
        next!(prog)
    end
    
    # compute stacking weights
    if label == "LSE"
        w = QP_stacking_weight(y_expect, y);
    else
        w = stacking_weight(lp_expect);
    end
    
    #return grid_all, weights, and label 
    return Dict(:grid_all => grid_all, :w => w, :label => label)
end

function sp_stacking_K_fold_MT(X, y, coords, deltasq_grid, phi_grid,
    nu_grid, Priors; K_fold = 10, seed = 123, label = "LSE", J = 300)
    ## experimential multithread version of sp_stacking_K_fold
    
    Random.seed!(seed);
    # pre-computation and pre-allocation #
    N = size(X, 2);
    CV_ind_ls = collect(Kfold(N, K_fold)); # index of train data in CV
    CV_ind_hold_ls = [setdiff(1:N, CV_ind_ls[k]) for k in 1:K_fold]; # index of held-out data in CV
    N_grid = length(deltasq_grid) * length(phi_grid) * length(nu_grid);
    nk_list = [length(x) for x in CV_ind_ls]; # be careful, different from the nk_list in my R code
    nk_k_list = [(N - x) for x in nk_list];   # This is the nk_list in my R code


    if X == Nothing()
        p = 0;
    else
        p = size(X, 1);
        Priors["inv_V_μ_β"] = Priors["inv_V_β"] * Priors["μβ"];
        XTX = X * X'; XTy = X * y;
        XTX_list = [XTX - X[:, CV_ind_hold_ls[k]] * X[:, CV_ind_hold_ls[k]]' for k in 1:K_fold];
        XTy_list = [XTy - X[:, CV_ind_hold_ls[k]] * y[CV_ind_hold_ls[k]] for k in 1:K_fold];
    end

    if label == "LSE"
        y_expect = Array{Float64, 2}(undef, N, N_grid);
    elseif label == "LP"
        lp_expect = Array{Float64, 2}(undef, N, N_grid);
        y_sq_sum_list = [norm(y[CV_ind_ls[k]])^2 for k in 1:K_fold];
    else 
        print("label has to be LSE or LP");
    end

    grid_phi_nu = vcat([[x y] for x in phi_grid, y in nu_grid]...);
    grid_all = vcat([[x y z] for x in phi_grid, y in nu_grid, z in deltasq_grid]...);
    L_grid_deltasq  = length(deltasq_grid);
    
    ## Compute expectation for stacking ##
    #println("start: ");
    #prog = Progress(size(grid_phi_nu, 1), 1, "Computing initial pass...", 50)
    Threads.@threads for i1 in 1:size(grid_phi_nu, 1)
        #println("Thread $(Threads.threadid()) handles $i1");
        #t_start = time()
        phi_pick = grid_phi_nu[i1, 1];
        nu_pick = grid_phi_nu[i1, 2];
        #Threads.@threads 
        for k in 1:K_fold
            Random.seed!(seed + (i1 - 1) * K_fold + k);
            if label == "LSE"
                y_expect[CV_ind_hold_ls[k], 
                    (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                stacking_prediction_LSE(coords, nu_pick, phi_pick, deltasq_grid, 
                    L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, 
                    nk_list, nk_k_list, y, X, XTX_list, XTy_list, Priors);
            else
                lp_expect[CV_ind_hold_ls[k], 
                    (i1 - 1) * L_grid_deltasq .+ (1:L_grid_deltasq)] = 
                stacking_prediction_LP(coords, nu_pick, phi_pick, deltasq_grid, 
                    L_grid_deltasq, k, CV_ind_ls, CV_ind_hold_ls, p, 
                    nk_list, nk_k_list, y, X, XTX_list, XTy_list,
                    y_sq_sum_list, Priors, J);
            end     
        end
        #println("Thread $(Threads.threadid()) finished $i1 in $(time() - t_start) seconds");
        #next!(prog)
    end
    
    # compute stacking weights
    if label == "LSE"
        w = QP_stacking_weight(y_expect, y);
    else
        w = stacking_weight(lp_expect);
    end
    
    #return grid_all, weights, and label 
    return Dict(:grid_all => grid_all, :w => w, :label => label)
end
=#

