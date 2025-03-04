module Statistics

export Histogram,
    bins_center

using Statistics

struct Histogram

    edges::AbstractRange{Float64}
    freqs::AbstractVector{UInt64}

    # Constructor from values
    function Histogram(values::AbstractVector{<:Real}, n_bins::Integer)
        low, high = extrema(values)
        bin_width = (high - low) / n_bins

        edges = range(low, high, length=n_bins + 1)

        freqs = zeros(UInt64, n_bins)
        for idx ∈ min.(floor.(UInt64, (values .- low) ./ bin_width) .+ 1, n_bins)
            freqs[idx] += 1
        end

        return new(edges, freqs)
    end

end

@inline bins_center(hist::Histogram) =
    hist.edges[begin:end-1] .+ (diff(hist.edges) ./ 2)

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

end
