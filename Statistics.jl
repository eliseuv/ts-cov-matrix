module Statistics

export Histogram

struct Histogram{T<:Real}
    edges::AbstractRange{T}
    freqs::AbstractVector{UInt64}

    # Constructor
    function Histogram(values::AbstractVector{<:Real}, n_bins::Integer)
        low, high = extrema(values)
        bin_width = (high - low) / n_bins

        edges = range(low, high, length=n_bins + 1)
        freqs = zeros(UInt64, n_bins)

        for idx ∈ min.(floor.(UInt64, (values .- low) ./ bin_width) .+ 1, n_bins)
            freqs[idx] += 1
        end

        return new{Float64}(edges, freqs)
    end

end

end
