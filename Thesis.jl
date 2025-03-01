module Thesis

export normalize_ts, normalize_ts!

using Statistics, LinearAlgebra

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
    cross_correlation_matrix(M::AbstractMatrix{<:Real})

Cross correlation matrix `G` of a given time series matrix `M_ts`.

    ``G = \frac{1}{N_{samples}} M_{ts}^T M_{ts}``

# Arguments:
- `M::AbstractMatrix`: `N×M` Matrix whose each of its `M` columns corresponds to a sample of a time series `Xₜ` of length `N`.
"""
@inline cross_correlation_matrix(M::AbstractMatrix{<:Real}) =
    Symmetric(M' * M) ./ size(M, 1)

end
