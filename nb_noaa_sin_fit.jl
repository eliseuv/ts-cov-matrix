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
    using Statistics, JLD2, Dates, DataFrames, Plots, LsqFit
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8da30af3-67dc-44f0-958e-57fdf85fa5e9
const df_ts = load_object("./data/noaa/1e5_uniform_points/sst_time_series.jld2")

# ╔═╡ 171f7d57-8c47-4d44-8b55-88cac95ba8ca
xdata = df_ts[!, :x5]

# ╔═╡ 0b905124-12cd-43bc-960e-3a8aa80b0059
tdata = 1:length(xdata)

# ╔═╡ 20998acf-e727-430a-a9b5-f1ec4d35569b
model(t, p) = p[1] .+ p[2] * sin.(((2π / 365) * t) .+ p[3])

# ╔═╡ 5ac714c3-7930-4ada-872c-cb8726d542f8
md"""
## Estimate initial values for the parameters
"""

# ╔═╡ 3350bb67-1a74-4cd4-b0f4-0030d7b02407
@doc"""
	approx_index(x::AbstractVector{<:Real}, x₀::Real)

Find index `i` for which `x[i] ≈ x₀`.

"""
@inline approx_index(x::AbstractVector{<:Real}, x₀::Real) = if x[begin] > x₀
	findfirst(<=(x₀), x)
else
	findfirst(>=(x₀), x)
end

# ╔═╡ 8a973528-b031-4237-947b-427aec9934b8
@doc"""
	init_params(x::AbstractVector{<:Real})

Estimate initial parameters for the model fit.
"""
function init_params(x::AbstractVector{<:Real})
	# Offset
	(xₘᵢₙ, xₘₐₓ) = extrema(x)
	B = (xₘᵢₙ + xₘₐₓ)  / 2

	# Amplitude
	A = (xₘₐₓ - xₘᵢₙ) / 2

	# Phase
	φ = approx_index(x, B) * 2π / 365
	
	return [B, A, φ]
end

# ╔═╡ 36cc3c24-8806-425a-bfed-08ee85cb4530
let x_init = model(tdata, init_params(xdata))
	plot(tdata, [xdata, x_init])
end

# ╔═╡ b138fa03-6375-4fcc-88cc-a9dc10f024f1
fit = curve_fit(model, tdata, xdata, init_params(xdata))

# ╔═╡ 9c133c00-aa7c-4db5-a8b0-0e44b22c3f9b
let xfit = model(tdata, fit.param)
	plot(tdata, [xdata, xfit])
end

# ╔═╡ f7019cbe-f124-463e-815c-345f9f401a83
function remove_seasonality(x::AbstractVector{<:Real})

	t = 1:length(x)

	fit = curve_fit(model, t, x, init_params(x))

	return x .- model(tdata, fit.param)
	
end

# ╔═╡ 815e0666-0122-4fcf-938e-58a52362dfb9
plot(tdata, [xdata .- mean(xdata), remove_seasonality(xdata)])

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8da30af3-67dc-44f0-958e-57fdf85fa5e9
# ╠═171f7d57-8c47-4d44-8b55-88cac95ba8ca
# ╠═0b905124-12cd-43bc-960e-3a8aa80b0059
# ╠═20998acf-e727-430a-a9b5-f1ec4d35569b
# ╟─5ac714c3-7930-4ada-872c-cb8726d542f8
# ╟─3350bb67-1a74-4cd4-b0f4-0030d7b02407
# ╠═8a973528-b031-4237-947b-427aec9934b8
# ╠═36cc3c24-8806-425a-bfed-08ee85cb4530
# ╠═b138fa03-6375-4fcc-88cc-a9dc10f024f1
# ╠═9c133c00-aa7c-4db5-a8b0-0e44b22c3f9b
# ╠═f7019cbe-f124-463e-815c-345f9f401a83
# ╠═815e0666-0122-4fcf-938e-58a52362dfb9
