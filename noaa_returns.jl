using DataFrames, Dates, GLM, JLD2

include("TimeSeries.jl")
include("HistStats.jl")

using .TimeSeries, .HistStats

@info "Loading dataframe"
const df_sst = let df = load("./data/noaa/random_points/uniform.jld2", "df_sst")
    t_min = minimum(Matrix(df[!, Not(:dates)]))
    # Offset temperature
    df[!, Not(:dates)] = df[!, Not(:dates)] .- t_min
    df
end
@info "Dataframe loaded"

const years_range = 1982:2024
const n_bins = 128
const x_min = 3e-2
const x_max = Inf
df = map(years_range) do y
    println(y)
    # Time series matrix
    M_ts = df_sst[year.(df_sst.dates).==y, Not(:dates)] |> Matrix{Float64}
    # Calculate returns
    rets = diff(M_ts, dims=1) ./ M_ts[begin:end-1, :] |>
           normalize_ts .|>
           abs |>
           filter(x -> isfinite(x) && x >= (1e-3))
    # Construct histogram
    hist = LogHistogram(rets, n_bins)
    # Fit
    df_fit = DataFrame(
        x=HistStats.log_midpoints(hist.edges),
        y=hist.freqs ./ sum(hist.freqs)
    ) |> filter(row -> row.y > 0)
    df_fit[!, :log_x] = log.(df_fit[!, :x])
    df_fit[!, :log_y] = log.(df_fit[!, :y])

    i_min = findfirst(>=(x_min), df_fit[!, :x])
    i_max = findlast(<=(x_max), df_fit[!, :x])

    ols = lm(@formula(log_y ~ log_x), df_fit[i_min:i_max, :])

    (intercept, slope) = coef(ols)

    (year=y, r2=r2(ols), intercept=intercept, slope=slope)

end |> DataFrame
@show df

save_object("./data/noaa/random_points/returns_plaw_fit.jld2", df)

