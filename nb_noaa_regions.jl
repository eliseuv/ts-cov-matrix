### A Pluto.jl notebook ###
# v0.20.21

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
    using Dates, DataFrames, JLD2, StatsBase, LinearAlgebra, GLM, CairoMakie, ColorSchemes
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 6a129f4d-77e5-42ee-a814-bc6bd1a2da7c
const years_range = 1982:2024

# ╔═╡ 34d791bc-332c-492d-8d10-47d230c3136d
const eigvals_low = Dict(load_object("./data/noaa/1e5_uniform_points/fft_filter/L=0/latitudes/lat_low_time_series_cov_eigvals.jld2"))

# ╔═╡ 5af3ac84-fe0e-4a19-bfa9-5d4539293487
const eigvals_high = Dict(load_object("./data/noaa/1e5_uniform_points/fft_filter/L=0/latitudes/lat_high_time_series_cov_eigvals.jld2"))

# ╔═╡ ca8e5d0a-9330-49b6-9b18-9544a9f4b821
eigvals_stats(eigvals) = let nbins = 128
	df = map(collect(eigvals)) do (y, eigvals)
		eigvals_hist = HistStats.Histogram(vec(eigvals), nbins)
		eigvals_mean = mean(eigvals_hist)
		eigvals_var = var(eigvals_hist)
		eigvals_max_mean = mean(eigvals[end, :])
		(year=y,
		 mean=eigvals_mean,
		 var=eigvals_var,
		 eigvals_max_mean=eigvals_max_mean)
	end |> DataFrame

	sort!(df, :year)
	
	df
end

# ╔═╡ 3cb9d0f2-ad66-4c1c-bff0-44d9025211ca
eigvals_low_stats = eigvals_stats(eigvals_low)

# ╔═╡ 64ec29f2-d0fa-48e0-98a4-5db4eba2a6fc
eigvals_high_stats = eigvals_stats(eigvals_high)

# ╔═╡ 8a7c2382-45f9-46b7-9605-56103d2a77d5
let fig = Figure(),
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (low latitudes)",
			  xlabel=L"n", ylabel=L"\langle\lambda_n\rangle",
			  #xscale=log10,
			  yscale=log10,
			  limits=((1, nothing), nothing)),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = mean(eigvals_low[year_sel], dims=2) |> vec |> reverse

		scatterlines!(ax, 1:length(λ), λ, color=c, markersize=3)
	
	end

	#axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	#save("./plots/noaa/spectra/eigval_hist_yearly.pdf", fig)
	
	fig
end

# ╔═╡ b76e06f6-b2ad-4e82-ba49-ab45724f3d9c
let fig = Figure(),
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (high latitudes)",
			  xlabel=L"n", ylabel=L"\langle\lambda_n\rangle",
			  #xscale=log10,
			  yscale=log10,
			  limits=((1, nothing), nothing)),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = mean(eigvals_high[year_sel], dims=2) |> vec |> reverse

		scatterlines!(ax, 1:length(λ), λ, color=c, markersize=3)
	
	end

	#axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	#save("./plots/noaa/spectra/eigval_hist_yearly.pdf", fig)
	
	fig
end

# ╔═╡ 5c8122f5-6c39-4da4-adfd-7dcdfc1da2c7
const (mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, 100)

# ╔═╡ b5d3e659-0134-4b55-8b61-bbd5be385d05
@inline clip_min(x_min, x) =
	if x < x_min
		x_min
	else
		x
	end

# ╔═╡ ba13ea9c-4649-49d3-854d-a2c785ef05cc
@inline clip_min(x_min) =
	x -> clip_min(x_min, x)

# ╔═╡ c5afafba-5ea3-4bd5-8a9f-e1d599f22a95
let fig = Figure(),
	eigvals = eigvals_low,
	years_range = years_range[begin]:6:years_range[end],
	λ_min = 0,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 1e-4,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (low latitudes)",
			  yscale=log10,
			  limits=((λ_min, λ_max),(y_min, 2))),
	nbins = 128,
	edges = range(λ_min, λ_max, length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	x = range(λ₋, λ_max, length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)
		#band!(ax, x, fill(y_min, length(y)), (clip_min(y_min).(y)); color = (c, 0.1))

		#I = sum(diff(hist.edges[1]) .* hist.weights)
		#println("$(year_sel): $(I)")
		
	end

	axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range)

	save("./plots/noaa/spectra/eigval_hist_yearly_low_lat.pdf", fig)
	
	fig
end

# ╔═╡ 760c1547-115f-4db7-a5bb-92480e0cb259
let fig = Figure(),
	eigvals = eigvals_high,
	years_range = years_range[begin]:6:years_range[end],
	λ_min = 0,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 1e-4,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (high latitudes)",
			  yscale=log10,
			  limits=((λ_min, λ_max),(y_min, 2))),
	nbins = 128,
	edges = range(λ_min, λ_max, length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	x = range(λ₋, λ_max, length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)
		#band!(ax, x, fill(y_min, length(y)), (clip_min(y_min).(y));color = (c, 0.1))

		#I = sum(diff(hist.edges[1]) .* hist.weights)
		#println("$(year_sel): $(I)")
		
	end

	axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range)

	save("./plots/noaa/spectra/eigval_hist_yearly_high_lat.pdf", fig)
	
	fig
end

# ╔═╡ de518cad-7789-48dc-9092-f7d47486f48c
let fig = Figure(),
	eigvals = eigvals_low,
	λ_min = 4e-4,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 3e-5,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (low latitudes)",
			  xscale=log10, yscale=log10,
			  limits=((λ_min, λ_max),(y_min, 2e2))),
	nbins = 128,
	edges = logrange(λ_min, λ_max, length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)
	
	x = logrange(λ₋, λ_max, length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		
		I = sum(diff(hist.edges[1]) .* y)
		#println("$(year_sel): $(I)")
	
	end

	axislegend(; merge = true)
	
	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	#save("./plots/noaa/spectra/eigval_hist_yearly_low_lat_loglog.pdf", fig)
	
	fig
end

# ╔═╡ fb04e65b-7cf5-4e10-9084-b4b2ce709df1
let fig = Figure(),
	eigvals = eigvals_high,
	λ_min = 4e-4,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 3e-5,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (high latitudes)",
			  xscale=log10, yscale=log10,
			  limits=((λ_min, λ_max),(y_min, 2e2))),
	nbins = 128,
	edges = logrange(λ_min, λ_max, length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)
	
	x = logrange(λ₋, λ_max, length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		
		I = sum(diff(hist.edges[1]) .* y)
		#println("$(year_sel): $(I)")
	
	end

	axislegend(; merge = true)
	
	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	#save("./plots/noaa/spectra/eigval_hist_yearly_high_lat_loglog.pdf", fig)
	
	fig
end

# ╔═╡ e4d32997-f322-4e7e-8d20-a5c8368b0fa6
let fig = Figure(),
	eigvals_ordered = reduce(hcat, map(years_range) do y
		map(mean, eachrow(eigvals_low[y]))
	end),
	ax = Axis(fig[1,1];
			  title="Eigenvalues Spatial Average",
			  ylabelsize=20,
			  ylabel=L"\langle \lambda_n \rangle", ylabelrotation=0,
			  yscale = log10,
			  limits=(extrema(years_range), nothing)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical = true)

	hlines!(ax, λ₋; color=:black, linewidth=2, label="Marchenko-Pastur Support")
	hlines!(ax, λ₊; color=:black, linewidth=2, label="Marchenko-Pastur Support")
	band!(ax, [extrema(years_range)...], fill(λ₋, 2), fill(λ₊, 2);
		  color=(:black, 0.3), label="Marchenko-Pastur Support")

	for (row, c) ∈ zip(eachrow(eigvals_ordered) |> reverse, cmap)

		scatterlines!(ax, years_range, row; color=c)

	end
	
	#scatterlines!(ax, eigvals_low_stats[!,:year], eigvals_low_stats[!, :var], color=:red)
	
	axislegend(; position=:lb, merge = true)
	
	save("./plots/noaa/spectra/eigval_ts_low_lat.pdf", fig)

	fig
end

# ╔═╡ ed59437f-d6ae-4fb6-8283-fd0574b0c259
let fig = Figure(),
	eigvals_ordered = reduce(hcat, map(years_range) do y
		map(mean, eachrow(eigvals_high[y]))
	end),
	ax = Axis(fig[1,1];
			  title="Eigenvalues Spatial Average",
			  ylabelsize=20,
			  ylabel=L"\langle \lambda_n \rangle", ylabelrotation=0,
			  yscale = log10,
			  limits=(extrema(years_range), nothing)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical = true)

	hlines!(ax, λ₋; color=:black, linewidth=2, label="Marchenko-Pastur Support")
	hlines!(ax, λ₊; color=:black, linewidth=2, label="Marchenko-Pastur Support")
	band!(ax, [extrema(years_range)...], fill(λ₋, 2), fill(λ₊, 2);
		  color=(:black, 0.3), label="Marchenko-Pastur Support")

	for (row, c) ∈ zip(eachrow(eigvals_ordered) |> reverse, cmap)

		scatterlines!(ax, years_range, row; color=c)

	end
	
	#scatterlines!(ax, eigvals_high_stats[!,:year], eigvals_high_stats[!, :var], color=:red)
	
	axislegend(; position=:lb, merge = true)
	
	save("./plots/noaa/spectra/eigval_ts_high_lat.pdf", fig)

	fig
end

# ╔═╡ c4ebaa95-0174-464d-9a5f-bbe6805b42df
 let fig = Figure(),
	eigvals = eigvals_low,
	λ_min = 1e-3,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 1e-3,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (low latitudes)",
			  yscale=log10,
			  xlabel=L"\lambda",
			  xtickformat = values -> [rich("10", superscript("$(round(Int, value))"))  for value in values],
			  limits=((log10(λ_min), log10(λ_max)),(y_min, 1e0))),
	nbins = 128,
	edges = range(log10(λ_min), log10(λ_max), length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	x = range(log10(λ₋), log(λ_max), length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(10 .^ x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(10 .^ x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")
	 
	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec .|> log10

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		
		I = sum(diff(hist.edges[1]) .* y)
		#println("$(year_sel): $(I)")
	
	end
	 
	axislegend(; merge = true)
	
	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	save("./plots/noaa/spectra/eigval_hist_yearly_low_lat_loglog.pdf", fig)
	
	fig
end

# ╔═╡ f04ec892-1640-4394-bc8e-f56ba2064958
 let fig = Figure(),
	eigvals = eigvals_high,
	λ_min = 1e-3,
	λ_max = maximum(y -> maximum(eigvals[y]), years_range),
	y_min = 1e-3,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum (high latitudes)",
			  yscale=log10,
			  xlabel=L"\lambda",
			  xtickformat = values -> [rich("10", superscript("$(round(Int, value))"))  for value in values],
			  limits=((log10(λ_min), log10(λ_max)),(y_min, 1e0))),
	nbins = 128,
	edges = range(log10(λ_min), log10(λ_max), length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	x = range(log10(λ₋), log(λ_max), length = 1000)
	lines!(ax, x, (clip_min(y_min) ∘ mp).(10 .^ x);
		   color=:black,
		   label="Marchenko-Pastur")
	band!(ax, x, fill(y_min, length(x)), (clip_min(y_min) ∘ mp).(10 .^ x);
		  color = (:black, 0.3),
		  label="Marchenko-Pastur")
	 
	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals[year_sel] |> filter(>=(0)) |> vec .|> log10

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		
		I = sum(diff(hist.edges[1]) .* y)
		#println("$(year_sel): $(I)")
	
	end

	axislegend(; merge = true)
	 
	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	save("./plots/noaa/spectra/eigval_hist_yearly_high_lat_loglog.pdf", fig)
	
	fig
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═6a129f4d-77e5-42ee-a814-bc6bd1a2da7c
# ╠═34d791bc-332c-492d-8d10-47d230c3136d
# ╠═5af3ac84-fe0e-4a19-bfa9-5d4539293487
# ╠═ca8e5d0a-9330-49b6-9b18-9544a9f4b821
# ╠═3cb9d0f2-ad66-4c1c-bff0-44d9025211ca
# ╠═64ec29f2-d0fa-48e0-98a4-5db4eba2a6fc
# ╠═8a7c2382-45f9-46b7-9605-56103d2a77d5
# ╠═b76e06f6-b2ad-4e82-ba49-ab45724f3d9c
# ╠═5c8122f5-6c39-4da4-adfd-7dcdfc1da2c7
# ╠═b5d3e659-0134-4b55-8b61-bbd5be385d05
# ╠═ba13ea9c-4649-49d3-854d-a2c785ef05cc
# ╠═c5afafba-5ea3-4bd5-8a9f-e1d599f22a95
# ╠═760c1547-115f-4db7-a5bb-92480e0cb259
# ╠═de518cad-7789-48dc-9092-f7d47486f48c
# ╠═fb04e65b-7cf5-4e10-9084-b4b2ce709df1
# ╠═e4d32997-f322-4e7e-8d20-a5c8368b0fa6
# ╠═ed59437f-d6ae-4fb6-8283-fd0574b0c259
# ╠═c4ebaa95-0174-464d-9a5f-bbe6805b42df
# ╠═f04ec892-1640-4394-bc8e-f56ba2064958
