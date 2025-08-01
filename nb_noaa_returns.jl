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
    using PlutoUI, Dates, JLD2, DataFrames, Statistics, StatsBase, GLM, CairoMakie
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

# ╔═╡ f47b0f8a-cfef-403b-82c5-d3d3eafe22ad
equinox_period(date::Date) =
	let y = year(date),
		eq1 = Date(y, 3, 20),
		eq2 = Date(y, 9, 23)
	
	if date < eq1
		(y-1, 2)
	elseif date < eq2
		(y, 1)
	else
		(y, 2)
	end
	
end

# ╔═╡ b7c5a862-46b1-4426-abb1-8810080bc394
const df_sst = let df = load(DATA_UNIFORM_PATH, "df_sst"),
	T_min = df[!, Not(:dates)] |> eachcol |> minimum |> minimum,
	T_kelvin = 273.15,
	ε = 1e-4
	# Offset by minimum temperature
	df[!, Not(:dates)] = df[!, Not(:dates)] .- T_min .+ ε
	# Label equionox periods
	insertcols!(df, 2, :equinox => equinox_period.(df[!, :dates]))
	df
end

# ╔═╡ b44f3d85-9c90-4e2c-b49a-065937107786
md"""
# Daily returns
"""

# ╔═╡ 6ac3583f-5a64-4d4a-8082-1d3d1cd6ae31
@inline calculate_returns(M_ts::AbstractMatrix{<:Real}) =
	diff(M_ts, dims=1) ./ M_ts[begin:end-1,:]

# ╔═╡ ab8f408e-3ccd-4828-8f66-756b7c65c1b6
md"""
# Select year
"""

# ╔═╡ f6f66349-ca6d-4dc8-b0c7-3e124b58b1ae
const years_range = 1982:2024 

# ╔═╡ 9f81977f-9cb5-4a1e-8b82-2a76d50f4b12
@bind year_sel Select(years_range)

# ╔═╡ f16807fb-0d9e-4e06-8d99-51b7ec378ef3
@inline get_time_series_matrix(y::Integer) =
	df_sst[year.(df_sst.dates).==y, r"x\d"] |> Matrix{Float64}

# ╔═╡ 5bf99978-14f0-4f19-81c3-52d016f43e77
const M_ts = get_time_series_matrix(year_sel)

# ╔═╡ c82dbf3c-69f8-46fe-b077-a94e9545d4fa
const daily_rets = calculate_returns(M_ts)

# ╔═╡ 6396d0f0-f109-434b-b8c5-d0031d7e1c49
let fig = Figure(),
	ax = Axis(fig[begin, begin]; yscale=log10)
	hist!(ax, daily_rets |> vec |> TimeSeries.normalize_ts, bins=128)

	fig
end

# ╔═╡ 4a3c0cff-6a3a-4078-92fd-fa9ba3f69480
md"""
# Seasons
"""

# ╔═╡ 92bee877-e2a4-4990-9bf1-000a94ef7c9f
md"""
# Period returns
"""

# ╔═╡ 64a69044-0bd8-4ff1-a60f-d296ee03c49a
function period_average(ts::AbstractVector, period::Integer)
	n = length(ts) ÷ period
	map(1:n) do k
		mean(@view ts[(1+(k-1)*period):k*period])
	end
end

# ╔═╡ d392f358-d5b6-4f53-8e09-8901bf8f8b01
period_average(M_ts::AbstractMatrix, period::Integer) =
	reduce(hcat, map(ts -> period_average(ts, period), eachcol(M_ts)))

# ╔═╡ 2c7d1def-f1d6-4c28-9d72-67d7baca0b0b
const period = 7

# ╔═╡ e2c9d3e2-3e8a-401e-9beb-dbd009b44709
const period_rets = period_average(M_ts, period) |> calculate_returns

# ╔═╡ ad8ab9e7-c544-42eb-a4a4-8d1aa2f8e5c6
let fig = Figure(),
	ax = Axis(fig[begin, begin]; yscale=log10)
	hist!(ax, period_rets |> TimeSeries.normalize_ts |> vec .|> abs, bins=128, normalization=:pdf)

	fig
end

# ╔═╡ fe1a95b3-b4ef-443a-a640-96fd232a5cda
@inline load_M_ts(y::Integer) =
	df_sst[year.(df_sst.dates).==y, Not(:dates)] |> Matrix{Float64}

# ╔═╡ e775cd4f-47b6-48e4-af78-a5d9bb82e7e4
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
	hist = HistStats.Histogram(rets, n_bins)
	#df_hist = DataFrame(
	#	x=HistStats.midpoints(hist.edges),
	#	y=hist.freqs ./ sum(hist.freqs)
	#) |> filter(row -> row.y > 0)

	#scatter!(ax, df_hist[!,:x], df_hist[!,:y])

	#save("./plots/returns/rets_$(y).pdf", fig)
	
end
	
end

# ╔═╡ Cell order:
# ╟─4ee071a6-3677-11f0-1088-f1da10eff286
# ╠═d7d2e436-3716-4b5b-8718-a7a2e827aa94
# ╟─958f37fb-8c3e-4233-83f2-994acc984227
# ╟─89216b8e-f977-4456-995f-5ccd0ec13447
# ╟─cd016e77-abf7-4076-8f02-f2757a69aebd
# ╟─8a54de12-fc66-4d84-bfef-7ea36d5efff2
# ╟─583a9c9f-aac8-47b4-934a-64bbdc921a9e
# ╟─f47b0f8a-cfef-403b-82c5-d3d3eafe22ad
# ╠═b7c5a862-46b1-4426-abb1-8810080bc394
# ╟─b44f3d85-9c90-4e2c-b49a-065937107786
# ╠═6ac3583f-5a64-4d4a-8082-1d3d1cd6ae31
# ╟─ab8f408e-3ccd-4828-8f66-756b7c65c1b6
# ╟─f6f66349-ca6d-4dc8-b0c7-3e124b58b1ae
# ╟─9f81977f-9cb5-4a1e-8b82-2a76d50f4b12
# ╠═f16807fb-0d9e-4e06-8d99-51b7ec378ef3
# ╠═5bf99978-14f0-4f19-81c3-52d016f43e77
# ╠═c82dbf3c-69f8-46fe-b077-a94e9545d4fa
# ╠═6396d0f0-f109-434b-b8c5-d0031d7e1c49
# ╟─4a3c0cff-6a3a-4078-92fd-fa9ba3f69480
# ╟─92bee877-e2a4-4990-9bf1-000a94ef7c9f
# ╠═64a69044-0bd8-4ff1-a60f-d296ee03c49a
# ╠═d392f358-d5b6-4f53-8e09-8901bf8f8b01
# ╠═2c7d1def-f1d6-4c28-9d72-67d7baca0b0b
# ╠═e2c9d3e2-3e8a-401e-9beb-dbd009b44709
# ╠═ad8ab9e7-c544-42eb-a4a4-8d1aa2f8e5c6
# ╠═fe1a95b3-b4ef-443a-a640-96fd232a5cda
# ╠═e775cd4f-47b6-48e4-af78-a5d9bb82e7e4
