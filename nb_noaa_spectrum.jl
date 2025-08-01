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
    using PlutoUI, Dates, DataFrames, JLD2, StatsBase, LinearAlgebra, GLM, CairoMakie, ColorSchemes
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 349dccc5-5805-4d3b-83cc-6bed00ebec4b
const eigvals_dict =  Dict(load_object("./data/noaa/1e5_uniform_points/sst_seasonality_removed_cov_eigvals_yearly.jld2"))

# ╔═╡ bda20425-1815-484f-8b54-4be5f6a828e3
const years_range = 1982:2024

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

# ╔═╡ caf1c8d4-8a45-4819-bffd-9911e65db467
let fig = Figure(),
	λ_min = 1e-6, λ_max = 2e2,
	ax = Axis(fig[1,1];
			  title="Covariance matrix spectrum",
			  xscale=log10, yscale=log10,
			  limits=((λ_min, λ_max),(1e-6, 1e4))),
	nbins = 128,
	edges = logrange(λ_min, λ_max, length=nbins+1),
	cbarPal = :jet,
	cmap = cgrad(colorschemes[cbarPal], length(years_range), categorical=true)

	for (year_sel, c) ∈ zip(years_range, cmap)

		λ = eigvals_dict[year_sel] |> filter(>=(0)) |> vec

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

	for (row, c) ∈ zip(eachrow(eigvals_ordered) |> reverse, cmap)

		scatterlines!(ax, years_range, row; color=c)

	end
	
	#scatterlines!(ax, df_stats[!,:year], df_stats[!, :mean], color=:black)
	
	#save("./plots/noaa/spectra/eigvals_ts.pdf", fig)

	fig
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═349dccc5-5805-4d3b-83cc-6bed00ebec4b
# ╠═bda20425-1815-484f-8b54-4be5f6a828e3
# ╠═8190ff39-0c10-47ac-81c7-d5597c18a9b3
# ╠═638f1abc-6d10-4444-8fcb-b325efd866e8
# ╠═31f0f0ee-d31c-4ccd-a015-d374de964483
# ╠═05427e18-91f6-4ab1-8dc2-93f5f3286646
# ╠═61fe4228-50a9-4eef-b1fe-827dd8a5e55b
# ╠═caf1c8d4-8a45-4819-bffd-9911e65db467
# ╠═e49aa19d-543a-442b-9f5f-c2be3bdbceb7
# ╠═27097599-64b3-41fa-968b-c424846649b6
