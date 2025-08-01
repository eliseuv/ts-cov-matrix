### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 545aaffc-f9e6-11ef-3071-a1ea565f5cb0
begin
    using Pkg
    Pkg.activate(".")
    include("/home/evf/.julia/pluto_notebooks/ingredients.jl")
end

# ╔═╡ 06f4e22b-9379-46a8-9477-824aa10f7bee
begin
    using Statistics, JLD2, Dates, DataFrames, Plots
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8da30af3-67dc-44f0-958e-57fdf85fa5e9
const df_ts = load_object("./data/noaa/1e5_uniform_points/sst_time_series.jld2")

# ╔═╡ 171f7d57-8c47-4d44-8b55-88cac95ba8ca
ts = df_ts[!, :x5]

# ╔═╡ eb0e5175-1cb8-4fe1-99a5-054606d0adcf
function get_seasonality(x::AbstractVector{<:Real}, period::Integer)

	n = length(x)

	map(1:period) do i
		mean(x[i:period:n])
	end
	
end

# ╔═╡ 63882865-f529-4c83-bcb2-4f8cf6d0cdf6
plot(1:365, get_seasonality(ts, 365))

# ╔═╡ d8c33cfe-c447-42df-b50b-39167ef7d16a
let period = 365,
	seasonality = get_seasonality(ts, period),
	n = length(ts),
	t = 1:n

	plot(t, [ts, repeat(seasonality, n ÷ period + 1)[begin:n]])
	
end

# ╔═╡ 725da8f4-85a6-40ef-b1f3-231e81bd0902
remove_seasonality(x::AbstractVector{<:Real}, period::Integer) =
	let n = length(x),
		t = 1:n

	x .- repeat(get_seasonality(x, period), n ÷ period + 1)[begin:n]
	
end

# ╔═╡ f4b59f48-9276-45f4-9022-09a10fd417a8
let period = 365,
	n = length(ts),
	t = 1:n

	plot(t, [ts, remove_seasonality(ts, period)])
	
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8da30af3-67dc-44f0-958e-57fdf85fa5e9
# ╠═171f7d57-8c47-4d44-8b55-88cac95ba8ca
# ╠═eb0e5175-1cb8-4fe1-99a5-054606d0adcf
# ╠═63882865-f529-4c83-bcb2-4f8cf6d0cdf6
# ╠═d8c33cfe-c447-42df-b50b-39167ef7d16a
# ╠═725da8f4-85a6-40ef-b1f3-231e81bd0902
# ╠═f4b59f48-9276-45f4-9022-09a10fd417a8
