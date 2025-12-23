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
    using PlutoUI, JLD2, DataFrames, Dates, LinearAlgebra, Statistics, StatsBase, Plots, CSV
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
end

# ╔═╡ 037d4d84-9900-4395-962a-a68a957d36af
const input_datafile = "./data/noaa/1e5_uniform_points/fft_filter/time_series_L=0.jld2"

# ╔═╡ d4ad916b-140f-413b-bb2e-97724793113f
df_sst = let df = load_object(input_datafile),
		     Tₘᵢₙ = df[!,Not(:dates)] |> (minimum ∘ (minimum .∘ eachrow)),
			 Tₖ = 273.15
	for x ∈ filter(n -> occursin(r"x\d+", n), names(df))
		df[!, x] = df[!, x] .+ Tₖ
	end
	df
end

# ╔═╡ 935fac35-dae3-4089-ade1-448bc0c43a3a
@inline get_time_series_matrix(y::Integer) =
	subset(df_sst, :dates => d -> year.(d) .== y)[!, Not(:dates)] |> Matrix

# ╔═╡ 812090cc-1489-48ce-be01-b25fc196054b
@inline calculate_returns(M_ts::AbstractMatrix{<:Real}) =
	(diff(M_ts, dims=1) ./ M_ts[begin:end-1,:]) |>
		vec |> filter(isfinite) |> TimeSeries.normalize_ts .|> abs

# ╔═╡ d0ba4675-8c78-4726-8a64-451824ecaf39
const years_range = 1982:2024

# ╔═╡ b51a10e4-0f9c-4297-abe1-0b822311c174
let years = 1982:6:2024,
	n = length(years),
	xₘᵢₙ= 1e-7, xₘₐₓ = 2e1,
	nbins = 512
	
	histogram(
		map(filter(<=(xₘₐₓ)) ∘ calculate_returns ∘ get_time_series_matrix, years),
		label=years .|> string |> permutedims,
		nbins=nbins,
		xscale=:log10, yscale=:log10,
		linewidth=2,
		line_z=(0:(n-1))', color=:viridis, colorbar=false,
		xlim=(xₘᵢₙ, xₘₐₓ)
	)

end

# ╔═╡ 3c709b49-0781-4721-8913-d5de8b4029b9
rets_dict = let years = 1982:6:2024,
	n = length(years),
	xₘᵢₙ = 1e-7, xₘₐₓ = 1e4,
	nbins = 128
	
	map(years) do y

		rets = get_time_series_matrix(y) |> calculate_returns |> filter(<=(xₘₐₓ))
		edges = logrange(xₘᵢₙ, xₘₐₓ, length=nbins+1)
		hist = fit(Histogram, rets, edges)
		
		y => (midpoints(edges), hist.weights)
	end

end

# ╔═╡ ecf5fa75-ebab-411e-8ab3-13ee65b39315
let plt = plot(; xscale=:log10)
	for (y, (x, w)) ∈ rets_dict
		(x, llw) = zip(filter(((_, w)) -> w > 0, zip(x, log.(w)) |> collect)...) |> collect
		#plot!(x, llw; label=y)
	end

	plt
end

# ╔═╡ ff83bf79-a4f4-466a-a4c4-9add5e042f49
for y ∈ years_range
	
	hist = fit(Histogram, get_time_series_matrix(y) |> calculate_returns |> filter(<=(20)); nbins=100000)

	df = DataFrame(bin_center=midpoints(hist.edges[begin]), counts=hist.weights)

	CSV.write("./rets_hist/rets_hist_y=$(y).csv", df)
	
	df
		
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═037d4d84-9900-4395-962a-a68a957d36af
# ╠═d4ad916b-140f-413b-bb2e-97724793113f
# ╠═935fac35-dae3-4089-ade1-448bc0c43a3a
# ╠═812090cc-1489-48ce-be01-b25fc196054b
# ╠═d0ba4675-8c78-4726-8a64-451824ecaf39
# ╠═b51a10e4-0f9c-4297-abe1-0b822311c174
# ╠═3c709b49-0781-4721-8913-d5de8b4029b9
# ╠═ecf5fa75-ebab-411e-8ab3-13ee65b39315
# ╠═ff83bf79-a4f4-466a-a4c4-9add5e042f49
