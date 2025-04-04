using JLD2, DataFrames, Dates, LinearAlgebra

include("TimeSeries.jl")
include("Statistics.jl")

using .TimeSeries, .Statistics

const df_sst = load("data/noaa/time_series/points.jld2", "df_sst")

get_matrix(name::AbstractString) =
    let df = dropmissing(select(df_sst, ["dates", name]), disallowmissing=true)
        reduce(hcat, sort(df[(Dates.year).(df.dates).==y, :], "dates")[begin:365, name] for y in 1982:2024)
    end

const points = names(df_sst) |> filter(!=("dates"))

const n_groups = 1_000;
const n_series = 100;
const n_bins = 256;

for point ∈ points[begin:begin]
    @show point
    M_ts = bootstrap_columns(get_matrix(point), 5, n_groups * n_series) |> normalize_ts
    @show size(M_ts)
    cov_eigvals = covariance_matrix_eigvals(M_ts, n_groups)
    @show size(cov_eigvals)
    hist = Histogram(vec(cov_eigvals), n_bins)
    @show hist

end
