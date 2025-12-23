module TimeSeries

export normalize_ts, normalize_ts!,
    covariance_matrix,
    covariance_matrix_eigvals,
    generate_columns,
    bootstrap_columns,
    split_columns,
    marchenko_pastur

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

"""
    persistence(x::AbstractVector{T}, ϵ::T=0) where {T<:Real}

Calculate the persistence of a time series `x`.

# Arguments
- `x::AbstractVector`: [TODO:description]
"""
@inline persistence(x::AbstractVector) =
    map(1:length(x)-1) do i
        if x[i+1] > x[i]
            findfirst(<=(x[i]), x[i+1:end])
        else
            findfirst(>=(x[i]), x[i+1:end])
        end
    end

@inline persistence(x::AbstractVector{T}, ϵ::T) where {T<:Real} =
    map(1:length(x)-1) do i
        if x[i+1] > x[i]
            findfirst(<=(x[i] - ϵ), x[i+1:end])
        else
            findfirst(>=(x[i] + ϵ), x[i+1:end])
        end
    end

@doc raw"""
    covariance_matrix(M::AbstractMatrix{<:Real})

Covariance matrix `G` of a given time series matrix `M`.

    ``G = \frac{1}{m} M^T M``
"""
@inline covariance_matrix(M::AbstractMatrix{<:Real}) =
    Symmetric(M' * M) ./ size(M, 1)

@doc raw"""
    covariance_matrix_eigvals(M::AbstractMatrix{<:Real})

Calculate the eigenvalues of the covariance matrix of a given time series matrix `M`.
"""
@inline covariance_matrix_eigvals(M::AbstractMatrix{<:Real}) =
    (eigvals ∘ covariance_matrix)(M)

@doc raw"""
    split_columns(M::AbstractMatrix, n_groups::Integer)

Split the columns of a given matrix `M` into `n_groups` groups resulting in a 3D array.
"""
@inline split_columns(M::AbstractMatrix, n_groups::Integer) =
    let (n_rows, n_cols) = size(M)
        @assert n_cols % n_groups == 0 "Number of columns $(n_cols) is not divisible by number of groups $(n_groups)"
        n_cols′ = n_cols ÷ n_groups
        reshape(M, (n_rows, n_cols′, n_groups))
    end

@doc raw"""
    covariance_matrix_eigvals(M::AbstractMatrix{<:Real}, n_groups::Integer)

Calculate the eigenvalues of the covariance matrix of a given time series matrix `M`.
"""
@inline covariance_matrix_eigvals(M::AbstractMatrix{<:Real}, n_groups::Integer) =
    dropdims(mapslices(covariance_matrix_eigvals, split_columns(M, n_groups), dims=(1, 2)), dims=(2,))

@doc raw"""
    generate_columns(M::AbstractMatrix{<:Real}, n_series::Integer, n_samples::Integer)

Generate `n_samples` columns of a matrix `M` by averaging random `n_series` columns.
"""
@inline generate_columns(M::AbstractMatrix{<:Real}, n_batch::Integer, n_samples::Integer) =
    let n_cols = size(M, 2)
        reduce(hcat, map(_ -> mean(M[:, randperm(n_cols)[begin:n_batch]], dims=2), 1:n_samples))
    end

@doc raw"""
    bootstrap_columns(M::AbstractMatrix{<:Real}, n_batch::Integer, n_samples::Integer)

Bootstrap sample the columns of a matrix `M` by averaging random `n_batch` columns `n_samples` times.
"""
@inline bootstrap_columns(M::AbstractMatrix{<:Real}, n_batch::Integer, n_samples::Integer) =
    let n_cols = size(M, 2)
        reduce(hcat, map(_ -> mean(M[:, randperm(n_cols)[begin:n_batch]], dims=2), 1:n_samples))
    end

@doc raw"""
    marchenko_pastur(mc_steps::Integer, n_samples::Integer)

Returns the Marchenko-Pastur distribution along with its support for the given time series matrix dimensions.
"""
@inline marchenko_pastur(mc_steps::Integer, n_samples::Integer) =
    let a = 1 + (n_samples / mc_steps),
        b = 2 * sqrt(n_samples / mc_steps),
        λ₋ = a - b,
        λ₊ = a + b

        f = function (λ)
            if λ₋ <= λ <= λ₊
                (mc_steps / (2 * π * n_samples)) * (sqrt((λ - λ₋) * (λ₊ - λ)) / λ)
            else
                0
            end
        end

        (f, (λ₋, λ₊))

    end

end
