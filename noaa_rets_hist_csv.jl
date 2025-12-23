using JLD2, DataFrames, Dates, LinearAlgebra, Statistics, StatsBase, CSV

include("TimeSeries.jl")
using .TimeSeries

const input_datafile = "./data/noaa/1e5_uniform_points/fft_filter/time_series_L=0.jld2"

df_sst = let df = load_object(input_datafile),
    Tₘᵢₙ = df[!, Not(:dates)] |> (minimum ∘ (minimum .∘ eachrow)),
    Tₖ = 273.15

    #= for x ∈ filter(n -> occursin(r"x\d+", n), names(df)) =#
    #=     df[!, x] = df[!, x] .+ Tₖ =#
    #= end =#

    df
end

show(df_sst)

@inline get_time_series_matrix(y::Integer) =
    subset(df_sst, :dates => d -> year.(d) .== y)[!, Not(:dates)] |> Matrix

@inline calculate_returns(M_ts::AbstractMatrix{<:Real}) =
    (diff(M_ts, dims=1) ./ M_ts[begin:end-1, :]) |>
    vec |> filter(isfinite) |> TimeSeries.normalize_ts .|> abs

const output_path = "./rets_hist_raw"
mkpath(output_path)

const years_range = 1982:2024
for y ∈ years_range

    hist = fit(Histogram, get_time_series_matrix(y) |> calculate_returns |> filter(<=(20)); nbins=100000)

    df = DataFrame(bin_center=midpoints(hist.edges[begin]), counts=hist.weights)
    show(df)

    CSV.write(joinpath(output_path, "rets_hist_y=$(y).csv"), df)

end
