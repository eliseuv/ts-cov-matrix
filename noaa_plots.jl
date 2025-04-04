using JLD2, DataFrames, Dates, CairoMakie, Printf

include("TimeSeries.jl")
include("HistStats.jl")

using .TimeSeries, .HistStats

# Plot theme
const theme = merge(theme_latexfonts(), Theme(fontsize=22))

# Load time series
const df_sst = load("data/noaa/time_series/points.jld2", "df_sst")
const years_range = 1982:2024

@inline get_matrix(name::AbstractString) =
    let df = dropmissing(select(df_sst, ["dates", name]), disallowmissing=true)
        reduce(hcat, sort(df[(Dates.year).(df.dates).==y, :], "dates")[begin:365, name] for y in years_range)
    end

const points = names(df_sst) |> filter(!=("dates"))

const n_groups = 1_000
const n_series = 100
const n_bins = 256

for point ∈ points
    @show point

    # Load time series matrix
    M_ts = get_matrix(point)

    # Bootstrapping
    M_ts_boot = TimeSeries.bootstrap_columns(M_ts, 5, n_groups * n_series) |> TimeSeries.normalize_ts

    # Calculate covariance matrix eigenvalues
    cov_eigvals = TimeSeries.covariance_matrix_eigvals(M_ts_boot, n_groups)

    # Fit histogram
    hist = HistStats.LogHistogram(filter(>(0), cov_eigvals), n_bins)

    with_theme(theme) do
        let fig = Figure(),
            λ = filter(>=(0), cov_eigvals) |> sort,
            hist = HistStats.LogHistogram(filter(>(0), λ), n_bins),
            x = HistStats.log_midpoints(hist.edges),
            ax = Axis(fig[1, 1];
                title="$(uppercasefirst(point)) Covariance matrix spectrum",
                xlabel=L"\lambda_i",
                xticks=1:50:300,
                xtickformat=indices -> [(@sprintf "%.1e" x[Int(idx)+1]) for idx ∈ indices],
                ylabel=L"\rho(\lambda)",
                limits=(nothing, (0, nothing)))

            barplot!(ax, hist.freqs, gap=0)

            save("./plots/$(point)_covariance_spectrum.pdf", fig)
        end
    end

end
