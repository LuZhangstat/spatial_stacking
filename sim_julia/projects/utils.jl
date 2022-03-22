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