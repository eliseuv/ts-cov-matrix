using JLD2, Statistics, DataFrames

include("Utils.jl")
include("HistStats.jl")
using .Utils, .HistStats

const input_datafile = "./data/noaa/1e5_uniform_points/fft_filter/time_series_L=1_cov_eigvals.jld2"
@show input_datafile

const output_datafile = append_to_filename(input_datafile, "_stats")
@show output_datafile

const nbins = 128
@show nbins

@info "Loading data..."
const cov_eigvals_pairs = load_object(input_datafile)

@info "Calculating stats..."
df = map(cov_eigvals_pairs) do (y, cov_eigvals)

    hist = Histogram(cov_eigvals, nbins)
    λ_mean = mean(hist)
    λ_var = var(hist, mean=λ_mean)

    λ_max = cov_eigvals[end, :]
    λ_max_mean = mean(λ_max)
    λ_max_var = varm(λ_max, λ_max_mean)

    (year=y,
        eigval_mean=λ_mean, eigval_var=λ_var,
        eigval_max_mean=λ_max_mean, eigval_max_var=λ_max_var)
end |> DataFrame
show(df)

@info "Saving..."
save_object(output_datafile, df)
