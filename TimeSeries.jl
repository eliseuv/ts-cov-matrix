module Thesis

export normalize_ts, normalize_ts!,
    covariance_matrix,
    covariance_matrix_eigvals,
    split_columns,
    bootstrap_columns

using Statistics, Random, LinearAlgebra

@doc raw"""
    normalize_ts(x::AbstractVector{<:Real})

Normalize a given time series vector `x`:

    ``xᵢ′ = (xᵢ - x̄) / sₓ``
"""
@inline function normalize_ts(x::AbstractVector{<:Real})
    x̄ = mean(x)
    sₓ = stdm(x, x̄)
    if sₓ != 0
        (x .- x̄) ./ sₓ
    else
        zeros(length(x))
    end
end

@doc raw"""
    normalize_ts!(x::AbstractVector{Float64})

Normalize a given time series vector `x` in place:

    ``xᵢ = (xᵢ - x̄) / sₓ``
"""
@inline function normalize_ts!(x::AbstractVector{Float64})
    x̄ = mean(x)
    sₓ = stdm(x, x̄)
    if sₓ != 0
        @views x .= (x .- x̄) ./ sₓ
    else
        fill!(x, 0)
    end
end

@doc raw"""
    normalize_ts(x::AbstractMatrix{<:Real})

Normalize each column `xᵢ` a given time series matrix `M`:

    ``xᵢ′ = (xᵢ - x̄) / sₓ``
"""
@inline normalize_ts(M::AbstractMatrix{<:Real}) =
    reduce(hcat, map(normalize_ts, eachcol(M)))


@doc raw"""
    normalize_ts!(x::AbstractMatrix{Float64})

Normalize in place each column `xᵢ` a given time series matrix `M`:

    ``xᵢ = (xᵢ - x̄) / sₓ``
"""
@inline normalize_ts!(M::AbstractMatrix{Float64}) =
    foreach(normalize_ts!, eachcol(M))

@doc raw"""
    covariance_matrix(M::AbstractMatrix{<:Real})

Covariance matrix `G` of a given time series matrix `M_ts`.

    ``G = \frac{1}{N_{samples}} M_{ts}^T M_{ts}``

# Arguments:
- `M::AbstractMatrix`: `N×M` Matrix whose each of its `M` columns corresponds to a sample of a time series `Xₜ` of length `N`.
"""
@inline covariance_matrix(M::AbstractMatrix{<:Real}) =
    Symmetric(M' * M) ./ size(M, 1)

@doc raw"""
    covariance_matrix_eigvals(M_ts::AbstractMatrix{<:Real})

Covariance matrix eigenvalues `G` of a given time series matrix `M_ts`.
"""
covariance_matrix_eigvals = eigvals ∘ covariance_matrix ∘ normalize_ts


@doc raw"""
    split_columns(M::AbstractMatrix, n_samples::Integer)

Split the columns of a given matrix `M` into `n_samples` groups resulting in a 3D array.
"""
@inline split_columns(M::AbstractMatrix, n_samples::Integer) =
    let (n_rows, n_cols) = size(M)
        @assert n_cols % n_samples == 0 "Number of columns must be divisible by number of samples"
        n_cols′ = n_cols ÷ n_samples
        reshape(M, (n_rows, n_cols′, n_samples))
    end

@doc raw"""
    bootstrap_columns(M::AbstractMatrix{<:Real}, n_batch::Integer, n_samples::Integer)

Bootstrap sample the columns of a matrix `M` by averaging random `n_batch` columns `n_samples` times.
"""
@inline bootstrap_columns(M::AbstractMatrix{<:Real}, n_batch::Integer, n_samples::Integer) =
    let n_cols = size(M, 2)
        reduce(hcat, map(_ -> mean(M[:, randperm(n_cols)[begin:n_batch]], dims=2), 1:n_samples))
    end

end
