### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 4ee071a6-3677-11f0-1088-f1da10eff286
begin
    using Pkg
    Pkg.activate(".")
    include("/home/evf/.julia/pluto_notebooks/ingredients.jl")
end

# ╔═╡ d7d2e436-3716-4b5b-8718-a7a2e827aa94
begin
    using PlutoUI, Dates, JLD2, DataFrames, Statistics, GLM, CairoMakie, StatsBase, CSV
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 958f37fb-8c3e-4233-83f2-994acc984227
md"""
# Load data
"""

# ╔═╡ 89216b8e-f977-4456-995f-5ccd0ec13447
const DATA_UNIFORM_PATH = "./data/noaa/random_points/uniform.jld2"

# ╔═╡ cd016e77-abf7-4076-8f02-f2757a69aebd
md"""
## Coordinates
"""

# ╔═╡ 8a54de12-fc66-4d84-bfef-7ea36d5efff2
const df_coords = load(DATA_UNIFORM_PATH, "df_coords")

# ╔═╡ 583a9c9f-aac8-47b4-934a-64bbdc921a9e
md"""
## Sea Surface Temperature
"""

# ╔═╡ b239324c-4aae-42fa-baca-8e0e8aa08e44
const df_sst = let df = load(DATA_UNIFORM_PATH, "df_sst"),
	t_min = minimum(Matrix(df[!, Not(:dates)]))
	# Offset temperature
	df[!, Not(:dates)] = df[!, Not(:dates)] .+ t_min
	df	
end

# ╔═╡ ab8f408e-3ccd-4828-8f66-756b7c65c1b6
md"""
# Select year
"""

# ╔═╡ f6f66349-ca6d-4dc8-b0c7-3e124b58b1ae
const years_range = 1982:2024 

# ╔═╡ 9f81977f-9cb5-4a1e-8b82-2a76d50f4b12
@bind year_sel Select(years_range)

# ╔═╡ f16807fb-0d9e-4e06-8d99-51b7ec378ef3
@inline load_M_ts(y::Integer) =
	df_sst[year.(df_sst.dates).==y, Not(:dates)] |> Matrix{Float64}

# ╔═╡ 5bf99978-14f0-4f19-81c3-52d016f43e77
const M_ts = df_sst[year.(df_sst.dates).==year_sel, Not(:dates)] |> Matrix{Float64}

# ╔═╡ d6534436-e27c-43e1-b70f-66e10d0a2a3f
md"""
# Returns analysis
"""

# ╔═╡ 5ca19c05-71c8-4873-a481-400f4c9eb931
@inline calculate_returns(ts::AbstractVector{<:Real}) = diff(ts) ./ ts[begin:end-1]

# ╔═╡ 6ac3583f-5a64-4d4a-8082-1d3d1cd6ae31
@inline calculate_returns(M_ts::AbstractMatrix{<:Real}) = diff(M_ts, dims=1) ./ M_ts[begin:end-1,:]

# ╔═╡ c82dbf3c-69f8-46fe-b077-a94e9545d4fa
const M_rets = calculate_returns(M_ts)

# ╔═╡ f48896ce-622d-4962-8226-160181e0f959
md"""
## Global Normalization
"""

# ╔═╡ ac36aa18-9991-4086-9954-07d22d8e6d49
let fig = Figure(),
	ax = Axis(fig[1,1];
			 yscale=log10),
	rets = calculate_returns(M_ts) |>
		filter(x -> isfinite(x)) |>
		TimeSeries.normalize_ts .|>
		abs
	
	hist!(ax, rets; bins=128, normalization=:pdf)

	fig
	
end

# ╔═╡ 2cbbf70b-7965-4201-98da-a5b7bbe5ee3b
md"""
## Individual Normalization
"""

# ╔═╡ 7d417b2d-54e8-4e49-971a-6e042c0674af
let fig = Figure(),
	ax = Axis(fig[1,1];
			 limits=(nothing, (0, 2.8))),
	rets = calculate_returns(M_ts) |>
		TimeSeries.normalize_ts .|>
		abs |>
		filter(x -> isfinite(x))

	hist = HistStats.Histogram(rets, 128)
	x = HistStats.midpoints(hist.edges)
	y = hist.freqs ./ sum(hist.freqs)
	
	loglogy = log.(-log.(y))
	logx = log.(x)
	
	scatter!(ax, logx, loglogy)

	fig
	
end

# ╔═╡ 43fe3e48-71a2-4c3f-9b1e-a9f715027c90
let fig = Figure(),
	ax = Axis(fig[1,1]; yscale=log10),
	k = 1,
	L = 365,
	σ = 0.1,
	n_series = 100_000
	
	x = range(0, 2π, length=L)
	
	M_ts = mapreduce(_ -> k .* sin.(x) .+ σ .* randn(L), hcat, 1:n_series)

	rets = calculate_returns(M_ts) |> TimeSeries.normalize_ts .|> abs |> vec
	
	hist!(ax, rets, bins=128)

	fig
	
end

# ╔═╡ 0b6d3ce9-b1a7-419e-9215-b1e9d7cc1c56


# ╔═╡ aae890dc-11fe-4360-93fb-893ee87ca7a9
const rets_ = calculate_returns(M_ts) |> TimeSeries.normalize_ts .|> abs |> vec

# ╔═╡ fff49a24-6387-41a8-ba5c-f7ec87c58b91
const rets = calculate_returns(M_ts) |>
	TimeSeries.normalize_ts .|>
	abs |>
    filter(x -> isfinite(x) && x >= (1e-3))

# ╔═╡ 1fba36a7-fc6e-4a75-8a82-f57aef291d2d
let fig = Figure(),
	ax = Axis(fig[1,1],
			  title = "Temperature Returns Histogram ($(year_sel))",
			  xlabel = L"| \Delta T |",
			  yscale=log10),
	n_bins = 128,
	x_min = 0.3,
	x_max = 5
	
	# Histogram
	hist = HistStats.Histogram(rets, n_bins)
	df_fit = DataFrame(
		x=HistStats.midpoints(hist.edges),
		y=hist.freqs ./ sum(hist.freqs)
	) |> filter(row -> row.y > 0)
	
	# Linear Fit
	df_fit[!, :log_x] = log.(df_fit[!, :x])
	df_fit[!, :log_y] = log.(df_fit[!, :y])
	i_min = findfirst(>=(x_min), df_fit[!,:x])
	i_max = findlast(<=(x_max), df_fit[!,:x])
	ols = lm(@formula(log_y ~ x), df_fit[i_min:i_max,:])
	(b, a) = coef(ols)
	rsq = round(r2(ols), digits=5)
	α =  round(a, digits=5)

	scatter!(ax, df_fit[!,:x], df_fit[!,:y])

	x = df_fit[i_min:i_max,:x]
	lines!(ax, x, exp(b) .* exp.(a .* x), color=:red)
	
	text!(ax , 1, 1, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"r^2 = %$(rsq)", fontsize=20)
	text!(ax , 1, 0.9, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"\alpha = %$(α)", fontsize=20)

	fig
	
end

# ╔═╡ 88714775-5a60-4283-826b-b930bcc6b046
let fig = Figure(),
	ax = Axis(fig[1,1],
			  title = "Temperature Returns Histogram ($(year_sel))",
			  xlabel = L"| \Delta T |")
	
	hist = fit(Histogram, rets)
	
end

# ╔═╡ c4150e96-c725-4cac-b833-9d9c21e6707d
let fig = Figure(),
	ax = Axis(fig[1,1],
			  title = "Temperature Returns Histogram ($(year_sel))",
			  xlabel = L"| \Delta T |", xscale=log10,
			  yscale=log10),
	n_bins = 128,
	x_min = 1e-2,
	x_max = 5e-1,
	# Histogram
	hist = HistStats.LogHistogram(rets, n_bins),
	df_hist = DataFrame(
		x=HistStats.log_midpoints(hist.edges),
		y=hist.freqs ./ sum(hist.freqs)
	) |> filter(row -> row.y > 0),
	# Linear Fit
	df_fit = DataFrame(
		log_x=log.(df_hist[!, :x]),
		log_y=log.(df_hist[!, :y])),
	i_min = findfirst(>=(x_min), df_hist[!,:x]),
	i_max = findlast(<=(x_max), df_hist[!,:x]),
	ols = lm(@formula(log_y ~ log_x), df_fit[i_min:i_max,:]),
	(b, a) = coef(ols),
	rsq = round(r2(ols), digits=5),
	α =  round(a, digits=5)

	scatter!(ax, df_hist[!,:x], df_hist[!,:y])

	x = df_hist[i_min:i_max,:x]
	lines!(ax, x, exp(b) .* (x .^ a), color=:red)
	
	text!(ax , 1, 1, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"r^2 = %$(rsq)", fontsize=20)
	text!(ax , 1, 0.9, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"\alpha = %$(α)", fontsize=20)

	fig
	
end

# ╔═╡ 81285255-699f-45b2-8898-68eaa58e63a0
# ╠═╡ disabled = true
#=╠═╡
for year_sel ∈ years_range
	
let fig = Figure(),
	ax = Axis(fig[1,1],
			  title = "Temperature Returns ($(year_sel))",
			  xlabel = L"| \Delta T |", xscale=log10,
			  yscale=log10),
	n_bins = 128,
	x_min = 2,
	x_max = Inf,
	rets = load_M_ts(year_sel) |>
		calculate_returns |>
	TimeSeries.normalize_ts .|>
	abs |>
    filter(x -> isfinite(x) && x >= (1e-3)),
	# Histogram
	hist = HistStats.LogHistogram(rets, n_bins),
	df_hist = DataFrame(
		x=HistStats.log_midpoints(hist.edges),
		y=hist.freqs ./ sum(hist.freqs)
	) |> filter(row -> row.y > 0),
	# Linear Fit
	df_fit = DataFrame(
		log_x=log.(df_hist[!, :x]),
		log_y=log.(df_hist[!, :y])),
	i_min = findfirst(>=(x_min), df_hist[!,:x]),
	i_max = findlast(<=(x_max), df_hist[!,:x]),
	ols = lm(@formula(log_y ~ log_x), df_fit[i_min:i_max,:]),
	(b, a) = coef(ols),
	rsq = round(r2(ols), digits=5),
	α =  round(a, digits=5)

	scatter!(ax, df_hist[!,:x], df_hist[!,:y])

	x = df_hist[i_min:i_max,:x]
	lines!(ax, x, exp(b) .* (x .^ a), color=:red)
	
	text!(ax , 1, 1, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"r^2 = %$(rsq)", fontsize=20)
	text!(ax , 1, 0.9, offset=(-4, -4), align = (:right, :top), space=:relative, text=L"\alpha = %$(α)", fontsize=20)

	save("./plots/noaa/returns/rets_$(year_sel).pdf", fig)
	
end

end
  ╠═╡ =#

# ╔═╡ 1d4182cc-0968-483f-a65c-f4a6ff840ef9
# ╠═╡ disabled = true
#=╠═╡
for y ∈ years_range

let fig = Figure(),
	ax = Axis(fig[1,1],
			  title = "Temperature Returns Histogram ($(y))",
			  xlabel = L"| \Delta T |",
			  yscale=log10,
			  limits=((0,1.5e-2), nothing)),
	rets = load_M_ts(y) |> calculate_returns .|> abs |> filter(>=(1e-3)),
	n_bins = 128,
	x_min = 2,
	x_max = 10,
	# Histogram
	hist = HistStats.Histogram(rets, n_bins),
	df_hist = DataFrame(
		x=HistStats.midpoints(hist.edges),
		y=hist.freqs ./ sum(hist.freqs)
	) |> filter(row -> row.y > 0)

	scatter!(ax, df_hist[!,:x], df_hist[!,:y])

	#save("./plots/returns/rets_$(y).pdf", fig)
	
end
	
end
  ╠═╡ =#

# ╔═╡ ad77b69e-6030-48cd-b868-36b7c25c4120
# ╠═╡ disabled = true
#=╠═╡
const df_ret_plaw = map(years_range) do y
	
	n_bins = 129
	x_min = 2
	x_max = Inf

	rets = load_M_ts(y) |> calculate_returns |> TimeSeries.normalize_ts .|> abs |> filter(>=(1e-3))

	# Histogram
	hist = HistStats.LogHistogram(rets, n_bins)
	df_hist = DataFrame(
		x=HistStats.log_midpoints(hist.edges),
		y=hist.freqs ./ sum(hist.freqs)
	) |> filter(row -> row.y > 0)
	# Linear Fit
	df_fit = DataFrame(
		log_x=log.(df_hist[!, :x]),
		log_y=log.(df_hist[!, :y]))
	i_min = findfirst(>=(x_min), df_hist[!,:x])
	i_max = findlast(<=(x_max), df_hist[!,:x])
	ols = lm(@formula(log_y ~ log_x), df_fit[i_min:i_max,:])
	(b, a) = coef(ols)

	(year=y, alpha=a, r2=r2(ols))
end |> DataFrame
  ╠═╡ =#

# ╔═╡ 75378692-7413-40a4-8dbd-b1d718050e92
#=╠═╡
let fig = Figure(),
	ax = Axis(fig[1,1],
			 title="Temperature Returns Power Law Exponent",
			  xticks=1982:6:2024,
			 limits=(nothing, (3, 4.5)))

	scatterlines!(ax, df_ret_plaw[!, :year], map(-, df_ret_plaw[!, :alpha]))
	
	save("./plots/noaa/returns/rets_plaw_exp.pdf", fig)
	
	fig
	
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─4ee071a6-3677-11f0-1088-f1da10eff286
# ╠═d7d2e436-3716-4b5b-8718-a7a2e827aa94
# ╟─958f37fb-8c3e-4233-83f2-994acc984227
# ╟─89216b8e-f977-4456-995f-5ccd0ec13447
# ╟─cd016e77-abf7-4076-8f02-f2757a69aebd
# ╟─8a54de12-fc66-4d84-bfef-7ea36d5efff2
# ╟─583a9c9f-aac8-47b4-934a-64bbdc921a9e
# ╠═b239324c-4aae-42fa-baca-8e0e8aa08e44
# ╟─ab8f408e-3ccd-4828-8f66-756b7c65c1b6
# ╟─f6f66349-ca6d-4dc8-b0c7-3e124b58b1ae
# ╟─9f81977f-9cb5-4a1e-8b82-2a76d50f4b12
# ╟─f16807fb-0d9e-4e06-8d99-51b7ec378ef3
# ╟─5bf99978-14f0-4f19-81c3-52d016f43e77
# ╟─d6534436-e27c-43e1-b70f-66e10d0a2a3f
# ╠═5ca19c05-71c8-4873-a481-400f4c9eb931
# ╠═6ac3583f-5a64-4d4a-8082-1d3d1cd6ae31
# ╠═c82dbf3c-69f8-46fe-b077-a94e9545d4fa
# ╟─f48896ce-622d-4962-8226-160181e0f959
# ╠═ac36aa18-9991-4086-9954-07d22d8e6d49
# ╟─2cbbf70b-7965-4201-98da-a5b7bbe5ee3b
# ╠═7d417b2d-54e8-4e49-971a-6e042c0674af
# ╠═43fe3e48-71a2-4c3f-9b1e-a9f715027c90
# ╠═0b6d3ce9-b1a7-419e-9215-b1e9d7cc1c56
# ╠═aae890dc-11fe-4360-93fb-893ee87ca7a9
# ╠═fff49a24-6387-41a8-ba5c-f7ec87c58b91
# ╠═1fba36a7-fc6e-4a75-8a82-f57aef291d2d
# ╠═88714775-5a60-4283-826b-b930bcc6b046
# ╠═c4150e96-c725-4cac-b833-9d9c21e6707d
# ╠═81285255-699f-45b2-8898-68eaa58e63a0
# ╠═1d4182cc-0968-483f-a65c-f4a6ff840ef9
# ╠═ad77b69e-6030-48cd-b868-36b7c25c4120
# ╠═75378692-7413-40a4-8dbd-b1d718050e92
