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

# ╔═╡ 545aaffc-f9e6-11ef-3071-a1ea565f5cb0
begin
    using Pkg
    Pkg.activate(".")
    include("/home/evf/.julia/pluto_notebooks/ingredients.jl")
end

# ╔═╡ 06f4e22b-9379-46a8-9477-824aa10f7bee
begin
    using PlutoUI, Statistics, Dates, StatsBase, LinearAlgebra, DataFrames, JLD2, CairoMakie, ColorSchemes
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ dc2197ef-b3cb-44e5-8b55-6735f14509f0
const cov_eigvals_dict = load_object("./data/noaa/1e5_uniform_points/fft_filter/time_series_L=0_cov_eigvals.jld2") |> Dict

# ╔═╡ a01a0b65-c55d-4463-89e0-de5dd55e6ec8
const years_range = keys(cov_eigvals_dict) |> collect |> sort

# ╔═╡ 0fea1672-e820-4532-99fb-349f5b0d93e5
const eigvals_ordered = reduce(hcat, map(sort(collect(cov_eigvals_dict), by=first)) do (_, ev)
	map(mean, eachrow(ev))
end)

# ╔═╡ bef73e62-a309-4205-9e13-a64db47ba401
let fig = Figure(),
	ax = Axis(fig[1,1];
			  title="Eigenvalues Spatial Average",
			  ylabelsize=20,
			  ylabel=L"\langle \lambda_n \rangle", ylabelrotation=0,
			  yscale = log10,
			  limits=(extrema(years_range), nothing)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical = true)

	for (row, c) ∈ zip(eachrow(eigvals_ordered) |> reverse, cmap)

		scatterlines!(ax, years_range, row; color=c)

	end
	
	#scatterlines!(ax, df_stats[!,:year], df_stats[!, :mean], color=:black)
	
	save("./plots/noaa/spectra/eigvals_ts.pdf", fig)

	fig
end

# ╔═╡ 145e76fb-efe8-4bbb-b9e3-e0cecfd90e27
const df_spectral_stats = map(years_range) do year_sel
	let eigvals = cov_eigvals_dict[year_sel],
		n_bins = 128,
		hist = HistStats.Histogram(eigvals[begin:end,:] |> vec, n_bins)
	λ_mean = mean(hist)
	λ_std = stdm(hist, λ_mean)
	λ_max = eigvals[end, :]
	λ_max_mean = mean(λ_max)
	λ_max_std = stdm(λ_max, λ_max_mean)
	λ_max_2_mean = mean(eigvals[end-1, :])
	(year=year_sel,
	 λ_mean=λ_mean, λ_std=λ_std,
	 λ_max_mean=λ_max_mean, λ_max_std=λ_max_std,
  	 λ_max_2_mean=λ_max_2_mean)
	end
end |> DataFrame

# ╔═╡ bf2fab74-4820-4b74-80ed-06c5db52703a
let fig = Figure(),
	ax = Axis(fig[1,1]; title="Average Eigenvalue",
			 xlabel="Year")

	scatterlines!(ax, df_spectral_stats[!,:year], df_spectral_stats[!, :λ_mean])

	save("./plots/noaa/spectra/eigval_avg.png", fig)
	
	fig
end

# ╔═╡ 0bbd8f42-6db9-45e5-8f47-53e899cb10ad
df_stats = load_object("./data/noaa/spatial_dist/global_stats.jld2")

# ╔═╡ 86ebb732-c637-4669-bf22-f21e1888bdf9
let fig = Figure(size=(600, 400))
	
	ax1 = Axis(fig[1,1];
			   yaxisposition=:right,
			   ylabel=L"\langle T \rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:red)
	
	ax2 = Axis(fig[1,1];
			   xlabel="Year",
			   ylabel=L"\langle \lambda \rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:blue)
	hidespines!(ax2)
	hidexdecorations!(ax2)
	
	
	scatterlines!(ax1, df_stats[!,:year], df_stats[!, :mean], color=:red)
	

	scatterlines!(ax2, df_spectral_stats[!,:year], df_spectral_stats[!, :λ_mean], color=:blue)
	
	save("./plots/noaa/spectra/eigval_avg.pdf", fig)

	fig
end

# ╔═╡ 8c0e2789-947e-433e-98bc-7f71ae8c316c
let fig = Figure(size=(600, 400))
	
	ax1 = Axis(fig[1,1];
			   yaxisposition=:right,
			   ylabel=L"\langle T \rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:red)
	
	ax2 = Axis(fig[1,1];
			   xlabel="Year",
			   ylabel=L"\langle \lambda_0 \rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:blue)
	hidespines!(ax2)
	hidexdecorations!(ax2)
	
	
	scatterlines!(ax1, df_stats[!,:year], df_stats[!, :mean], color=:red)
	

	scatterlines!(ax2, df_spectral_stats[!,:year], df_spectral_stats[!, :λ_max_mean], color=:blue)
	
	save("./plots/noaa/spectra/eigval_max_avg.pdf", fig)

	fig
end

# ╔═╡ de8e2638-c5ed-459f-a8cd-6fb420b8f7c9
@bind year_sel Select(years_range)

# ╔═╡ 9d21ece1-9f1f-4724-883e-1e3e8a57aed1
const cov_eigvals = cov_eigvals_dict[year_sel]

# ╔═╡ 2dfc8e72-feda-4079-b409-05a24460967e
let fig = Figure(),
	ax = Axis(fig[1,1],
			 yscale=log10),
	λ = cov_eigvals[begin:end, :] |> vec

	hist!(ax, λ;
		  bins=128,
		  normalization=:pdf)

	(mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, 100)
	x = range(λ₋, λ₊, 1000)
	lines!(ax, x, mp.(x), color=:red)

	fig
end

# ╔═╡ 54fdd604-f165-411d-9631-01a7351bc06c
let fig = Figure(),
	ax = Axis(fig[1,1],
			 yscale=log10)
	λ = cov_eigvals[begin:end, :] |> vec

	hist!(ax, λ;
		  bins=128,
		  normalization=:pdf)

	(mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, 100)
	x = range(λ₋, λ₊, 1000)
	lines!(ax, x, mp.(x), color=:red)

	fig
end

# ╔═╡ 8c4600fe-a05a-434a-8c97-61e60390b472
let fig = Figure(),
	ax = Axis(fig[1,1],
				xscale=log10, yscale=log10,
			 limits=((1e-7, 1e5),(1e-8, 1e4)))
	λ = cov_eigvals[begin:end, :] |> filter(>=(0)) |> vec
	nbins = 128

	hist = fit(Histogram, λ, logrange(extrema(λ)..., length=nbins+1)) |> normalize

	x = midpoints(hist.edges[1])
	y = hist.weights

	scatterlines!(ax, x, y)

	fig
end

# ╔═╡ 79b8d713-52a3-4083-a593-700655cb4695
for year_sel ∈ years_range
	
	fig = Figure()
	ax = Axis(fig[1,1],
				xscale=log10, yscale=log10,
			 limits=((1e-7, 1e5),(1e-8, 1e4)))
	nbins = 128
	
	λ = cov_eigvals_dict[year_sel] |> filter(>=(0)) |> vec

	hist = fit(Histogram, λ, logrange(extrema(λ)..., length=nbins+1)) |> normalize

	x = midpoints(hist.edges[1])
	y = hist.weights

	scatterlines!(ax, x, y)
	
	#(mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, 100)
	#x = range(λ₋, λ₊, 1000)
	#lines!(ax, x, mp.(x))

	save("./plots/noaa/spectra/eigvals_hist_$(year_sel).pdf", fig)
	
end

# ╔═╡ 789c012e-8877-4838-ab5c-dc262955afc0
let fig = Figure(),
	λ_min = 1e-7, λ_max = 1e5,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum",
			  xscale=log10, yscale=log10,
			  limits=((λ_min, λ_max),(1e-8, 1e4))),
	nbins = 128,
	edges = logrange(λ_min, λ_max, length=nbins+1),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = cov_eigvals_dict[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c, label = "$(year_sel)")
	
	end

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	save("./plots/noaa/spectra/eigval_hist_yearly.pdf", fig)
	
	fig
end

# ╔═╡ Cell order:
# ╟─545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═dc2197ef-b3cb-44e5-8b55-6735f14509f0
# ╟─a01a0b65-c55d-4463-89e0-de5dd55e6ec8
# ╠═0fea1672-e820-4532-99fb-349f5b0d93e5
# ╠═bef73e62-a309-4205-9e13-a64db47ba401
# ╟─145e76fb-efe8-4bbb-b9e3-e0cecfd90e27
# ╠═bf2fab74-4820-4b74-80ed-06c5db52703a
# ╠═0bbd8f42-6db9-45e5-8f47-53e899cb10ad
# ╠═86ebb732-c637-4669-bf22-f21e1888bdf9
# ╠═8c0e2789-947e-433e-98bc-7f71ae8c316c
# ╟─de8e2638-c5ed-459f-a8cd-6fb420b8f7c9
# ╠═9d21ece1-9f1f-4724-883e-1e3e8a57aed1
# ╠═2dfc8e72-feda-4079-b409-05a24460967e
# ╠═54fdd604-f165-411d-9631-01a7351bc06c
# ╠═8c4600fe-a05a-434a-8c97-61e60390b472
# ╠═79b8d713-52a3-4083-a593-700655cb4695
# ╠═789c012e-8877-4838-ab5c-dc262955afc0
