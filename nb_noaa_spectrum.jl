### A Pluto.jl notebook ###
# v0.20.18

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
    using PlutoUI, Dates, DataFrames, JLD2, StatsBase, LinearAlgebra, GLM, CairoMakie, ColorSchemes, CSV
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ fc57f6fd-9586-42ad-9f93-ea36f0d9dc58
@inline clip_min(x_min, x) =
	if x < x_min
		x_min
	else
		x
	end

# ╔═╡ bd306009-0f72-48c1-9079-8adb874e86ae
@inline clip_min(x_min) =
	x -> clip_min(x_min, x)

# ╔═╡ b0db9423-0e58-4a42-8bab-3b0f49357250
const (mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, 100)

# ╔═╡ 4322a95f-4718-43ed-886a-37684dddec79
const input_datafile = "./data/noaa/1e5_uniform_points/fft_filter/L=0/time_series_L=0_cov_eigvals.jld2"

# ╔═╡ 349dccc5-5805-4d3b-83cc-6bed00ebec4b
const eigvals_dict =  Dict(load_object(input_datafile))

# ╔═╡ bda20425-1815-484f-8b54-4be5f6a828e3
const years_range = 1982:2024

# ╔═╡ bf9db7a4-f870-432d-8605-ec8d3b6f0593
let fig = Figure(),
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum",
			  xlabel=L"n", ylabel=L"\langle\lambda_n\rangle",
			  xscale=log10,
			  yscale=log10,
			  limits=((1, nothing), nothing)),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = mean(eigvals_dict[year_sel], dims=2) |> vec |> reverse

		scatterlines!(ax, 1:length(λ), λ, color=c, markersize=3)
	
	end

	#axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	#save("./plots/noaa/spectra/eigval_hist_yearly.pdf", fig)
	
	fig
end

# ╔═╡ 0b5662ca-c2d6-4765-9630-6816e99c16fa
let output_dir = "./send/sst_eigvals/",
	nbins = 256

	mkpath(output_dir)

	for year_sel ∈ years_range
		println(year_sel)
		
		λ = eigvals_dict[year_sel] |> filter(>=(0)) |> vec

		edges = range(extrema(λ)..., nbins + 1)
		hist = fit(Histogram, λ, edges) |> normalize

		df = DataFrame(lambda=midpoints(hist.edges[1]),
					  dist=hist.weights)

		I = sum(diff(hist.edges[1]) .* hist.weights)
		println("$(year_sel): $(I)")

		CSV.write(joinpath(output_dir, "$(year_sel).csv"), df)
	end

end

# ╔═╡ 6f46378c-5dd7-4581-a0e8-76fe6d081595
let fig = Figure(),
	years_range = years_range[begin]:6:years_range[end],
	λ_min = 0,
	λ_max = maximum(y -> maximum(eigvals_dict[y]), years_range),
	y_min = 3e-5,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum",
			  yscale=log10,
			  limits=((λ_min, λ_max),(y_min, 2))),
	nbins = 100,
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

		λ = eigvals_dict[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		I = sum(diff(hist.edges[1]) .* hist.weights)
		println("$(year_sel): $(I)")
		
	end

	axislegend(; merge = true)

	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	save("./plots/noaa/spectra/eigval_hist_yearly.pdf", fig)
	
	fig
end

# ╔═╡ caf1c8d4-8a45-4819-bffd-9911e65db467
let fig = Figure(),
	λ_min = 4e-4,
	λ_max = maximum(y -> maximum(eigvals_dict[y]), years_range),
	y_min = 3e-5,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum",
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

		λ = eigvals_dict[year_sel] |> filter(>=(0)) |> vec

		hist = fit(Histogram, λ, edges) |> normalize

		x = midpoints(hist.edges[1])
		y = hist.weights

		lines!(ax, x, y, color=c)

		
		I = sum(diff(hist.edges[1]) .* y)
		println("$(year_sel): $(I)")
	
	end

	axislegend(; merge = true)
	
	Colorbar(fig[begin, end+1], colormap=cmap,
			 limits=extrema(years_range), nsteps=length(years_range), 
			 ticks=years_range[begin]:6:years_range[end])

	save("./plots/noaa/spectra/eigval_hist_yearly_loglog.pdf", fig)
	
	fig
end


# ╔═╡ e0d2e1e6-0a90-4eee-844d-3e880fb052bd
	const df_stats = load_object("./data/noaa/1e5_uniform_points/fft_filter/time_series_L=0_cov_eigvals_stats.jld2")

# ╔═╡ 2f150511-d25a-4d49-ab91-03ff967c67e3
let fig = Figure(),
	ax = Axis(fig[1,1],
			  title="Average eigenvalue",
			  xlabel="year",
			  ylabel=L"\langle\lambda\rangle", ylabelrotation=0,
			  limits=(extrema(years_range), nothing))

	scatterlines!(ax, df_stats[!,:year], df_stats[!,:eigval_mean])

	save("./plots/noaa/spectra/eigval_avg.pdf", fig)
	
	fig
	
end

# ╔═╡ 548b08c1-b27e-49b2-887b-d5bb80169f33
df_temp_stats = load_object("./data/noaa/spatial_dist/global_stats.jld2")

# ╔═╡ af344fe6-8db5-4d98-ab9b-3423b667807d
let fig = Figure(size=(600, 400)),
	
	ax1 = Axis(fig[1,1],
			   title="Average largest eigenvalue",
			   xlabel="Year",
			   yaxisposition=:right,
			   ylabel=L"\langle T \rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:red,
			   limits=(extrema(years_range), nothing))

	scatterlines!(ax1, df_temp_stats[!,:year], df_temp_stats[!,:mean], color=:red)
	
	ax2 = Axis(fig[1,1],
			   yaxisposition=:left,
			   ylabel=L"\langle\lambda_\text{max}\rangle",
			   ylabelsize=20, ylabelrotation=0,
			   ylabelcolor=:blue,
			   limits=(extrema(years_range), nothing))

	scatterlines!(ax2, df_stats[!,:year], df_stats[!,:eigval_max_mean], color=:blue)
	
	save("./plots/noaa/spectra/eigval_max_avg.pdf", fig)
	
	fig
	
end

# ╔═╡ b41d28b2-ae0c-405d-9a2b-17c67ce3e18f


# ╔═╡ 8190ff39-0c10-47ac-81c7-d5597c18a9b3
@bind period_sel Select(years_range)

# ╔═╡ 638f1abc-6d10-4444-8fcb-b325efd866e8
cov_eigvals = eigvals_dict[period_sel]

# ╔═╡ 31f0f0ee-d31c-4ccd-a015-d374de964483
let fig = Figure(),
	ax = Axis(fig[1,1];
			  yscale = log10)

	hist!(ax, cov_eigvals |> vec; bins=128)

	fig
		  
end

# ╔═╡ 05427e18-91f6-4ab1-8dc2-93f5f3286646
let fig = Figure(),
	λ = cov_eigvals |> vec |> filter(>(0)) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i", xscale=log10,
        ylabel=L"\lambda_i", yscale=log10)

    scatter!(ax, λ)

    fig
end

# ╔═╡ 61fe4228-50a9-4eef-b1fe-827dd8a5e55b
let fig = Figure(),
	λ = cov_eigvals |> filter(>=(0)) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance matrix spectrum CDF",
        xlabel=L"\lambda_i", xscale=log10,
        ylabel=L"F(\lambda)",
		limits=(extrema(λ), nothing))

    scatterlines!(ax, λ, range(0, 1, length(λ)))

    fig
end

# ╔═╡ e49aa19d-543a-442b-9f5f-c2be3bdbceb7
const eigvals_ordered = reduce(hcat, map(years_range) do y
	map(mean, eachrow(eigvals_dict[y]))
end)

# ╔═╡ 27097599-64b3-41fa-968b-c424846649b6
let fig = Figure(),
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
	
	#scatterlines!(ax, df_stats[!,:year], df_stats[!, :mean], color=:black)
	
	axislegend(; position=:lb, merge = true)
	
	save("./plots/noaa/spectra/eigvals_ts.pdf", fig)

	fig
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═fc57f6fd-9586-42ad-9f93-ea36f0d9dc58
# ╠═bd306009-0f72-48c1-9079-8adb874e86ae
# ╠═b0db9423-0e58-4a42-8bab-3b0f49357250
# ╠═4322a95f-4718-43ed-886a-37684dddec79
# ╠═349dccc5-5805-4d3b-83cc-6bed00ebec4b
# ╠═bda20425-1815-484f-8b54-4be5f6a828e3
# ╠═bf9db7a4-f870-432d-8605-ec8d3b6f0593
# ╠═0b5662ca-c2d6-4765-9630-6816e99c16fa
# ╠═6f46378c-5dd7-4581-a0e8-76fe6d081595
# ╠═caf1c8d4-8a45-4819-bffd-9911e65db467
# ╠═27097599-64b3-41fa-968b-c424846649b6
# ╠═e0d2e1e6-0a90-4eee-844d-3e880fb052bd
# ╠═2f150511-d25a-4d49-ab91-03ff967c67e3
# ╠═548b08c1-b27e-49b2-887b-d5bb80169f33
# ╠═af344fe6-8db5-4d98-ab9b-3423b667807d
# ╠═b41d28b2-ae0c-405d-9a2b-17c67ce3e18f
# ╠═8190ff39-0c10-47ac-81c7-d5597c18a9b3
# ╟─638f1abc-6d10-4444-8fcb-b325efd866e8
# ╠═31f0f0ee-d31c-4ccd-a015-d374de964483
# ╠═05427e18-91f6-4ab1-8dc2-93f5f3286646
# ╟─61fe4228-50a9-4eef-b1fe-827dd8a5e55b
# ╠═e49aa19d-543a-442b-9f5f-c2be3bdbceb7
