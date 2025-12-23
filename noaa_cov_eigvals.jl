@doc"""
    This scripts takes the SST time series and calculates the covariance matrix eigenvalues for each year.
"""

using JLD2, DataFrames, Dates

include("Utils.jl")
include("TimeSeries.jl")
include("HistStats.jl")

using .Utils, .TimeSeries, .HistStats

@inline equinox_year(date::Date) =
    let y = year(date),
        equinox_date = Date(y, 9, 23)

        if date < equinox_date
            y
        else
            y + 1
        end

    end

const input_datafile = "./data/noaa/1e5_uniform_points/fft_filter/L=0/latitudes/lat_high_time_series.jld2"
@show input_datafile

const output_datafile = append_to_filename(input_datafile, "_cov_eigvals")
@info output_datafile

@info "Loading data..."
df = load_object(input_datafile)
show(df)

df[!, :year] = year.(df[!, :dates])

const n_total = count(s -> occursin(r"x\d+", s), names(df))
@show n_total

const n_series = 100
@show n_series

const n_groups = n_total ÷ n_series
@show n_groups

const n_cols = n_groups * n_series
@show n_cols

@info "Calculating eigenvalues..."
const cov_eigvals = [
    key.year =>
        TimeSeries.covariance_matrix_eigvals(Matrix{Float64}(gdf[!, r"x\d+"])[:, 1:n_cols] |>
                                             normalize_ts,
            n_groups)
    for (key, gdf) ∈ pairs(groupby(df, :year))
]

save_object(output_datafile, cov_eigvals)
