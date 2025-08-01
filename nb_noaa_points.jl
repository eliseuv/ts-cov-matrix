### A Pluto.jl notebook ###
# v0.20.6

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
	using PlutoUI, JLD2, DataFrames, Dates, CairoMakie, ColorSchemes, Printf
	TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
	HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 89c874a6-0d42-4199-8fff-49e3896ecd4b
const df_sst = load("data/noaa/time_series/points.jld2", "df_sst")

# ╔═╡ 5b79ee54-9a17-4013-886a-101ea7c03aa4
const years_range = 1982:2024

# ╔═╡ 609bbf46-9414-476d-b8bc-bc374fc5c543
@inline get_matrix(name::AbstractString) =
    let df = dropmissing(select(df_sst, ["dates", name]), disallowmissing=true)
        reduce(hcat, sort(df[(Dates.year).(df.dates).==y, :], "dates")[begin:365, name] for y in years_range)
	end

# ╔═╡ 579b98fd-232b-4b05-b8b0-6981dabde405
const points = names(df_sst) |> filter(!=("dates"))

# ╔═╡ 8030e3bf-386a-4fb2-a4ce-114e3e3d6b44
const n_groups = 1_000

# ╔═╡ d7094df0-753c-4cfa-84c6-1fd1e5e77932
const n_series = 100

# ╔═╡ ece27cdf-321b-4655-bcd2-d5bbd2832fbf
const n_bins = 256

# ╔═╡ f7665fd1-69a4-481e-acbe-06e4cf7769ce
@bind point Select(points)

# ╔═╡ ab42151e-4791-472a-b9f4-e41c6e0990d3
md"""
# $(uppercasefirst(point))
"""

# ╔═╡ 2f2092f3-ef11-42e5-aede-e4ec1a3c9f95
md"""
# Time Series
"""

# ╔═╡ f0307466-be5d-4aca-8b4e-8533607ece22
const period = 1:365

# ╔═╡ e5396610-6f12-4dc7-b7f8-809a40ff6afb
M_ts = get_matrix(point)[period,:]

# ╔═╡ a29b8c8c-4c22-474b-b76c-93ede23604e3
let fig = Figure(),
    ax = Axis(fig[1,1];
        title="Time Series",
        xlabel="day",
        ylabel="Temperature (°C)",
        limits=((0, length(period)), nothing)),
    cbarPal = :viridis

    cmap = cgrad(colorschemes[cbarPal], size(M_ts,2), categorical = true)
    for (ts, color) in zip(eachcol(M_ts), cmap)
        lines!(ax, ts, color=color)
    end
    Colorbar(fig[begin, end+1], limits=extrema(years_range), colormap=cmap, ticks=years_range[begin]:3:years_range[end])

    fig
end

# ╔═╡ 27274e6a-f044-4f30-aa1c-466754c8ef4b
md"""
# Return analysis
"""

# ╔═╡ 6f655eeb-e951-44ce-a101-9e72ad834e87
let fig = Figure(),
	ax = Axis(fig[1,1],
			  yscale=Makie.pseudolog10),
	M_ts_ret = diff(M_ts, dims=1) ./ M_ts[begin:end-1,:]
	
	hist!(ax,filter(>=(0), vec(M_ts_ret)), bins=200)

	fig
end

# ╔═╡ 7bdd3df7-b78b-4cc8-94c5-cca6ba2f64db
let M_ts_ret = diff(M_ts, dims=1) ./ M_ts[begin:end-1,:],
	rets = filter(>=(0), vec(M_ts_ret)),
	hist = HistStats.Histogram(rets, 200)

	x = HistStats.midpoints(hist.edges)
	y = hist.freqs
	
end

# ╔═╡ 489526fa-6e3f-4a65-9526-ddb39c4f89bb
let fig = Figure(),
	ax = Axis(fig[1,1],
			  yscale=Makie.pseudolog10),
	M_ts_ret = diff(M_ts, dims=1) ./ M_ts[begin:end-1,:]
	
	hist!(ax, abs.(vec(M_ts_ret)), bins=200)

	fig
end

# ╔═╡ 5f897937-ab00-4c18-981e-2672e792876f
md"""
# Time Series Bootstrapping
"""

# ╔═╡ 7df00aff-4b09-4876-ba92-3651ca3b2cde
M_ts_boot = TimeSeries.bootstrap_columns(M_ts, 5, n_groups * n_series) |> TimeSeries.normalize_ts

# ╔═╡ 1df1dcc4-0a63-49a2-8612-59351a36e141
let fig = Figure(),
	nsamples=100,
    ax = Axis(fig[1,1];
        title="Normalized bootstrapped Time Series ($(nsamples) samples)",
        xlabel="day",
        ylabel="Temperature (°C)",
        limits=((0, length(period)), nothing)),
    cbarPal = :viridis,
	M = @view M_ts_boot[:,begin:nsamples]

    cmap = cgrad(colorschemes[cbarPal], size(M,2), categorical = true)
    for (ts, color) in zip(eachcol(M), cmap)
        lines!(ax, ts, color=color)
    end

    fig
end

# ╔═╡ b59c18c4-ecdb-4c61-8fa5-52d62a42e173
cov_eigvals = TimeSeries.covariance_matrix_eigvals(M_ts_boot, n_groups)

# ╔═╡ 1c31e891-c38f-4088-b5a4-d0f3fe4675f4
cov_eigvals[end,:]

# ╔═╡ 6b67b02e-957a-46a1-80c4-1b844189d52d
let fig = Figure(),
	λ = cov_eigvals |> vec |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i",
        ylabel=L"\lambda_i")

    scatter!(ax, λ)

    fig
end

# ╔═╡ 9afdfc55-fba2-4881-8da9-f8a911c1ed5d
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i",
        ylabel=L"\lambda_i", yscale=log10)

    scatter!(ax, λ)

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

	fig
end

# ╔═╡ a131aeeb-11cd-4036-a348-05af2b1e9969
function marchenko_pastur((λ₋, λ₊), mc_steps, n_samples)
	return function(λ)
		(mc_steps/(2*π*n_samples)) * (sqrt((λ - λ₋) * (λ₊ - λ)) / λ)
	end
end

# ╔═╡ d88b065f-1d27-4b82-9419-35d5d0ff1fb3
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
	hist = HistStats.Histogram(filter(>(0), λ), n_bins*50),
	x = HistStats.midpoints(hist.edges),
	mp = marchenko_pastur(extrema(λ), length(period), n_series)
	ax = Axis(fig[1,1];
        title="Covariance matrix spectrum",
        xlabel=L"\lambda_i",
        ylabel=L"\rho(\lambda)", yscale=Makie.pseudolog10,
		limits=((0, 100), (0, nothing)),
		)

    barplot!(ax, x, hist.freqs, gap=0)
	lines!(ax, x, mp.(x))

	

	fig
end

# ╔═╡ Cell order:
# ╟─545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╟─06f4e22b-9379-46a8-9477-824aa10f7bee
# ╟─89c874a6-0d42-4199-8fff-49e3896ecd4b
# ╟─5b79ee54-9a17-4013-886a-101ea7c03aa4
# ╠═609bbf46-9414-476d-b8bc-bc374fc5c543
# ╟─579b98fd-232b-4b05-b8b0-6981dabde405
# ╟─8030e3bf-386a-4fb2-a4ce-114e3e3d6b44
# ╟─d7094df0-753c-4cfa-84c6-1fd1e5e77932
# ╟─ece27cdf-321b-4655-bcd2-d5bbd2832fbf
# ╟─f7665fd1-69a4-481e-acbe-06e4cf7769ce
# ╟─ab42151e-4791-472a-b9f4-e41c6e0990d3
# ╟─2f2092f3-ef11-42e5-aede-e4ec1a3c9f95
# ╠═f0307466-be5d-4aca-8b4e-8533607ece22
# ╠═e5396610-6f12-4dc7-b7f8-809a40ff6afb
# ╠═a29b8c8c-4c22-474b-b76c-93ede23604e3
# ╟─27274e6a-f044-4f30-aa1c-466754c8ef4b
# ╠═6f655eeb-e951-44ce-a101-9e72ad834e87
# ╠═7bdd3df7-b78b-4cc8-94c5-cca6ba2f64db
# ╠═489526fa-6e3f-4a65-9526-ddb39c4f89bb
# ╟─5f897937-ab00-4c18-981e-2672e792876f
# ╟─7df00aff-4b09-4876-ba92-3651ca3b2cde
# ╟─1df1dcc4-0a63-49a2-8612-59351a36e141
# ╟─b59c18c4-ecdb-4c61-8fa5-52d62a42e173
# ╟─1c31e891-c38f-4088-b5a4-d0f3fe4675f4
# ╟─6b67b02e-957a-46a1-80c4-1b844189d52d
# ╟─9afdfc55-fba2-4881-8da9-f8a911c1ed5d
# ╟─2b5addbf-31b2-4ebc-aa03-6de53bdc0c07
# ╟─8f5404e0-df0b-4c82-9aea-1dbc41faba0e
# ╠═d88b065f-1d27-4b82-9419-35d5d0ff1fb3
# ╠═a131aeeb-11cd-4036-a348-05af2b1e9969
