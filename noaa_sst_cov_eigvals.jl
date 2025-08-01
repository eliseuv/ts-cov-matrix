@doc"""
    This scripts takes the SST time series and calculates the covariance matrix eigenvalues for each year.
"""

using JLD2, DataFrames, Dates

include("TimeSeries.jl")
include("HistStats.jl")

using .TimeSeries, .HistStats

@inline equinox_year(date::Date) =
    let y = year(date),
        equinox_date = Date(y, 9, 23)

        if date < equinox_date
            y
        else
            y + 1
        end

    end

df = load_object("./data/noaa/1e5_uniform_points/sst_time_series_seasonality_removed.jld2")

df[!, :year] = year.(df[!, :dates])

const n_total = count(s -> occursin(r"x\d+", s), names(df))
const n_series = 100
@assert n_total % n_series == 0
const n_groups = n_total ÷ n_series

@show n_total n_series n_groups

const cov_eigvals = [
    key.year => TimeSeries.covariance_matrix_eigvals(Matrix{Float64}(gdf[!, r"x\d+"]), n_groups)
    for (key, gdf) ∈ pairs(groupby(df, :year))
]

save_object("./data/noaa/1e5_uniform_points/sst_seasonality_removed_cov_eigvals_yearly.jld2", cov_eigvals)
