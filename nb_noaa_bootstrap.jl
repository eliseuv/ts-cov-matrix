### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 545aaffc-f9e6-11ef-3071-a1ea565f5cb0
begin
	using Pkg
	Pkg.activate(".")
	include("/home/evf/.julia/pluto_notebooks/ingredients.jl")
end

# ╔═╡ 06f4e22b-9379-46a8-9477-824aa10f7bee
begin
	using PlutoUI, Statistics, JLD2, DataFrames, Dates, CairoMakie, ColorSchemes, Printf, GLM
	TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
	HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 90b1261a-2868-42e7-bbd6-9c089687db8c
@inline celsius_to_kelvin(T::Real) = T + 273.15

# ╔═╡ 07421c76-f4e5-45d7-b56b-dfcf1738c16f
md"""
# Map
"""

# ╔═╡ a0cd91d0-f500-4874-a3af-c84925c60fb3
const DATA_PATH = "./data/noaa/random_points/uniform.jld2"

# ╔═╡ fe2a37a5-5f3d-4601-aa28-9af20df5938c
const df_coords = load(DATA_PATH, "df_coords")

# ╔═╡ 41498e2c-1285-49f7-a031-986e525c4151
const df_sst = load(DATA_PATH, "df_sst")

# ╔═╡ 386cdb94-ac6d-43e0-bf87-b2c69397c84a
let fig = Figure(),
	temp_avg = reduce(+, eachcol(df_sst[!, Not(:dates)])) ./ ncol(df_sst),
    ax = Axis(fig[1,1];
        title="Average Temperature",
        xlabel="day",
        ylabel="Temperature (°C)")

	lines!(ax, df_sst[!, :dates], temp_avg)

	save("./plots/avg_temperature.pdf", fig)
	
	fig
	
end

# ╔═╡ 2f2092f3-ef11-42e5-aede-e4ec1a3c9f95
md"""
# Time Series
"""

# ╔═╡ 5b79ee54-9a17-4013-886a-101ea7c03aa4
const years_range = 1982:2024

# ╔═╡ 2362e3d0-df4d-4f1d-9448-89b6ecbee7ce
@bind year Select(years_range)

# ╔═╡ c8ea2780-8f5f-49f2-a3ed-d81a5e40033f
df_year = df_sst[(Dates.year).(df_sst.dates) .== year, :]

# ╔═╡ 17418b05-0fb4-49ce-8380-3223ed163363
let fig = Figure(),
	dates = df_year[!, :dates],
    ax = Axis(fig[1,1];
        title="Time Series",
        xlabel="day",
        ylabel="Temperature (°C)")
	
    for ts in eachcol(select(df_year, Not(:dates)))[begin:32]
        lines!(ax, dates, ts)
    end

	save("./plots/ts_example_y$(year).pdf", fig)

    fig
end

# ╔═╡ cd9f4ebb-8f8f-4f2a-9858-a9eaa69279bb
const M_ts = convert(Array{Float64}, reduce(hcat, eachcol(select(df_year, Not(:dates)))))

# ╔═╡ 5f897937-ab00-4c18-981e-2672e792876f
md"""
# Single Year Time Series Bootstrapping
"""

# ╔═╡ 8030e3bf-386a-4fb2-a4ce-114e3e3d6b44
const n_groups = 1_000

# ╔═╡ d7094df0-753c-4cfa-84c6-1fd1e5e77932
const n_series = 100

# ╔═╡ f4dfcd4d-35e7-4ec2-955f-7b515bf9e2e8
const M_ts_boot = TimeSeries.bootstrap_columns(M_ts, 5, n_groups*n_series) |> TimeSeries.normalize_ts

# ╔═╡ 1df1dcc4-0a63-49a2-8612-59351a36e141
let fig = Figure(),
	nsamples=100,
    ax = Axis(fig[1,1];
        title="Normalized bootstrapped Time Series ($(nsamples) samples)",
        xlabel="day",
        ylabel="Temperature (°C)",
        limits=((0, size(M_ts_boot, 1)), nothing)),
    cbarPal = :viridis,
	M = @view M_ts_boot[:,begin:nsamples]

    cmap = cgrad(colorschemes[cbarPal], size(M,2), categorical = true)
    for (ts, color) in zip(eachcol(M), cmap)
        lines!(ax, ts, color=color)
    end

	save("./plots/ts_bootstrapped_y$(year).pdf", fig)
	
    fig
end

# ╔═╡ b59c18c4-ecdb-4c61-8fa5-52d62a42e173
cov_eigvals = TimeSeries.covariance_matrix_eigvals(M_ts_boot, n_groups)

# ╔═╡ 6b67b02e-957a-46a1-80c4-1b844189d52d
let fig = Figure(),
	λ = cov_eigvals |> vec |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i",
        ylabel=L"\lambda_i")

    scatter!(ax, λ)

	
	save("./plots/zipf_y$(year).pdf", fig)

    fig
end

# ╔═╡ 9afdfc55-fba2-4881-8da9-f8a911c1ed5d
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i", xscale=log10,
        ylabel=L"\lambda_i")

	logλ = log.(λ)
    scatter!(ax, log.(logλ .- minimum(logλ) .+ 0.1))
	
	save("./plots/zipf_y$(year)_log.pdf", fig)

    fig
end

# ╔═╡ 2b5addbf-31b2-4ebc-aa03-6de53bdc0c07
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance matrix spectrum CDF",
        xlabel=L"\lambda_i", xscale=log10,
        ylabel=L"F(\lambda)",
		limits=(extrema(λ), (0, nothing)))

    scatterlines!(ax, λ, range(0, 1, length(λ)))

    fig
end

# ╔═╡ ece27cdf-321b-4655-bcd2-d5bbd2832fbf
const n_bins = 256

# ╔═╡ 8f5404e0-df0b-4c82-9aea-1dbc41faba0e
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
	hist = HistStats.LogHistogram(filter(>(0), λ), n_bins),
	x = HistStats.log_midpoints(hist.edges),
	ax = Axis(fig[1,1];
        title="Covariance matrix spectrum",
        xlabel=L"\log(\lambda_i)",
		xticks=1:50:300,
		xtickformat= indices -> [ (@sprintf "%.2e" x[Int(idx)+1]) for idx ∈ indices],
        ylabel=L"\rho(\lambda)", yscale=Makie.pseudolog10,
		limits=(nothing, (0, nothing)),
		)

    barplot!(ax, hist.freqs, gap=0)

	save("./plots/eigvals_hist_y$(year).pdf", fig)

	fig
end

# ╔═╡ d88b065f-1d27-4b82-9419-35d5d0ff1fb3
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
	(mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, n_series)
	ax = Axis(fig[1,1];
        title="Covariance matrix spectrum",
        xlabel=L"\lambda_i",
        ylabel=L"\rho(\lambda)", yscale=Makie.log10,
		limits=((0, 60), (0.001, nothing)),
		)

    #barplot!(ax, x, hist.freqs, gap=0)
	hist!(ax, λ, bins=256, normalization=:pdf)

	x = 0:0.00001:5
	lines!(ax, x, mp.(x), color=:red)
	
	save("./plots/eigvals_hist_y$(year)_mp.pdf", fig)

	fig
end

# ╔═╡ 50285faf-7fbd-4f6e-98d4-d36b477e079d
# ╠═╡ disabled = true
#=╠═╡
@inline get_matrix(year::Integer) = 
	 reduce(hcat, eachcol(select(df_sst[(Dates.year).(df_sst.dates) .== year, :], Not(:dates))))	
  ╠═╡ =#

# ╔═╡ c9e8172d-ffa7-4fe7-8866-cc30c314dbcb
#=╠═╡
cov_eigvals_dict = map(years_range) do year
	M_ts = get_matrix(year)
	M_ts_boot = TimeSeries.bootstrap_columns(M_ts, 10, n_groups*n_series) |> TimeSeries.normalize_ts
	cov_eigvals = TimeSeries.covariance_matrix_eigvals(M_ts_boot, n_groups)
	year => cov_eigvals
end
  ╠═╡ =#

# ╔═╡ 31ced615-f68a-4ea4-b88f-c845d69c88aa
#=╠═╡
df_year_stats = map(cov_eigvals_dict) do (year, cov_eigvals)
	hist = HistStats.Histogram(cov_eigvals, n_bins)
	(year=year, mean=mean(hist), var=var(hist), eigval_max=mean(cov_eigvals[end,:]))
end |> DataFrame
  ╠═╡ =#

# ╔═╡ 641710b5-aa38-406e-9895-f0c19afb4e97
#=╠═╡
let fig = Figure(),
	ax = Axis(fig[1,1], 
			 title="Yearly mean eigenvalue bootstrapped from 1000 uniformly distributed points",
			 xlabel="year")

	scatterlines!(ax, df_year_stats[!,:year], df_year_stats[!, :var])

	#save("./plots/yearly_mean_eigvals.pdf", fig)
	
	fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═90b1261a-2868-42e7-bbd6-9c089687db8c
# ╟─07421c76-f4e5-45d7-b56b-dfcf1738c16f
# ╟─a0cd91d0-f500-4874-a3af-c84925c60fb3
# ╟─fe2a37a5-5f3d-4601-aa28-9af20df5938c
# ╟─41498e2c-1285-49f7-a031-986e525c4151
# ╟─386cdb94-ac6d-43e0-bf87-b2c69397c84a
# ╟─2f2092f3-ef11-42e5-aede-e4ec1a3c9f95
# ╟─5b79ee54-9a17-4013-886a-101ea7c03aa4
# ╟─2362e3d0-df4d-4f1d-9448-89b6ecbee7ce
# ╟─c8ea2780-8f5f-49f2-a3ed-d81a5e40033f
# ╠═17418b05-0fb4-49ce-8380-3223ed163363
# ╟─cd9f4ebb-8f8f-4f2a-9858-a9eaa69279bb
# ╟─5f897937-ab00-4c18-981e-2672e792876f
# ╟─8030e3bf-386a-4fb2-a4ce-114e3e3d6b44
# ╟─d7094df0-753c-4cfa-84c6-1fd1e5e77932
# ╟─f4dfcd4d-35e7-4ec2-955f-7b515bf9e2e8
# ╟─1df1dcc4-0a63-49a2-8612-59351a36e141
# ╟─b59c18c4-ecdb-4c61-8fa5-52d62a42e173
# ╠═6b67b02e-957a-46a1-80c4-1b844189d52d
# ╠═9afdfc55-fba2-4881-8da9-f8a911c1ed5d
# ╠═2b5addbf-31b2-4ebc-aa03-6de53bdc0c07
# ╟─ece27cdf-321b-4655-bcd2-d5bbd2832fbf
# ╠═8f5404e0-df0b-4c82-9aea-1dbc41faba0e
# ╠═d88b065f-1d27-4b82-9419-35d5d0ff1fb3
# ╠═50285faf-7fbd-4f6e-98d4-d36b477e079d
# ╠═c9e8172d-ffa7-4fe7-8866-cc30c314dbcb
# ╠═31ced615-f68a-4ea4-b88f-c845d69c88aa
# ╠═641710b5-aa38-406e-9895-f0c19afb4e97
