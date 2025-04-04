module HistStats

export
    Histogram, LogHistogram,
    midpoints, log_midpoints

using Statistics

struct Histogram{T<:Real}

    edges::AbstractVector{T}
    freqs::AbstractVector{UInt64}

    function Histogram(vals::AbstractArray{<:Real}, edges::AbstractVector{T}) where {T<:Real}
        freqs = zeros(UInt64, length(edges) - 1)
        for idx ∈ map(x -> findfirst(edges[begin:end-1] .<= x .< edges[begin+1:end]), vals) |> filter(!isnothing)
            freqs[idx] += 1
        end

        return new{T}(edges, freqs)
    end

    function Histogram(vals::AbstractArray{<:Real}, n_bins::Integer)
        low, high = extrema(vals)
        bin_width = (high - low) / n_bins

        edges = range(low, high, length=n_bins + 1)

        freqs = zeros(UInt64, n_bins)
        for idx ∈ min.(floor.(UInt64, (vals .- low) ./ bin_width) .+ 1, n_bins)
            freqs[idx] += 1
        end

        return new{Float64}(edges, freqs)
    end

end

@inline LogHistogram(vals::AbstractArray{<:Real}, n_bins::Integer) =
    let (low, high) = extrema(vals)
        Histogram(vals, logrange(low, high, n_bins + 1))
    end

@inline Statistics.mean(hist::Histogram) =
    let n = sum(hist.freqs)
        sum(e * (f / n) for (e, f) ∈ zip(hist.edges, hist.freqs))
    end

@inline function Statistics.var(hist::Histogram; corrected::Bool=true, mean=nothing)
    mean_used = if mean === nothing
        Statistics.mean(hist)
    elseif isa(mean, Real)
        mean
    else
        throw(ArgumentError("invalid value of mean, $(mean)::$(typeof(mean))"))
    end

    n = sum(hist.freqs)
    return sum(abs2(e - mean_used) * (f / (n - Int(corrected))) for (e, f) ∈ zip(hist.edges, hist.freqs))

end

@inline midpoints(x::AbstractVector{<:Real}) =
    @views (x[begin:end-1] .+ x[begin+1:end]) ./ 2

@inline log_midpoints(x::AbstractVector{<:Real}) =
    @views sqrt.(x[begin:end-1] .* x[begin+1:end])

end
