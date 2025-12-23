using JLD2, DataFrames, Dates, FFTW

@inline approx_idx(ran, val) =
    round(Int, (val - first(ran)) / step(ran) + 1)

@inline function band_stop_filter!(F, ω, ω₀, L=0)
    i = approx_idx(ω, ω₀)
    F[i-L:i+L] .= 0
end

@inline function apply_filter(x)
    # Number of samples
    Nₛ = length(x)

    # Sampling rate
    # Daily => Average number of days per year
    fₛ = 365.25

    # Frequencies
    ω = rfftfreq(Nₛ, fₛ)


    # Fourier transform
    F = rfft(x)

    nₘₐₓ = floor(Int, ω[end])
    for n ∈ 1:nₘₐₓ
        band_stop_filter!(F, ω, n, 1)
    end

    # Inverse Fourier transform
    x′ = irfft(F, Nₛ)

    return x′

end


const input_datafile = "./data/noaa/1e5_uniform_points/time_series.jld2"
@show input_datafile
const output_datafile = "./data/noaa/1e5_uniform_points/fft_filter/time_series_L=1.jld2"
@show output_datafile

@info "Loading data..."
df = filter(:dates => d -> year(d) ∈ 1982:2024,
    load_object(input_datafile))
show(df)

@info "Applying filter..."
mapcols!(apply_filter ∘ Vector{Float64}, df, cols=Not(:dates))
show(df)

@info "Saving..."
save_object(output_datafile, df)
