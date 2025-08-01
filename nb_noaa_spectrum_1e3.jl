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
	using PlutoUI, Random, Statistics, JLD2, DataFrames, Dates, CairoMakie, ColorSchemes, Printf, GLM
	TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
	HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8f814715-3830-42be-85d4-1f29c2b0bffb
const DATA_PATH = "./data/noaa/random_points/uniform.jld2"

# ╔═╡ b3dc1a6f-769e-43b3-ac16-f2d23e8cb71f
const df_coords = load(DATA_PATH, "df_coords")

# ╔═╡ 8002063f-c2a7-40d8-b8c8-8c244ecec2e9
const df_sst = let df = load(DATA_PATH, "df_sst"),
				   t_min = minimum(Matrix(df[!, Not(:dates)]))

	df[!, Not(:dates)] = df[!, Not(:dates)] .- t_min
	df
	
end

# ╔═╡ 3f0226c8-bce6-4c30-992c-6a2c61cd5dd1
const years_range = 1982:2024

# ╔═╡ 70c494fb-b3f8-432d-8948-ff1e690aeaa2
@bind y Select(years_range)

# ╔═╡ 9b4958f0-a21e-47a5-a246-9e05fe952f90
df_year = df_sst[(Dates.year).(df_sst.dates) .== y, :]

# ╔═╡ b8212742-d884-4819-8f32-4fe9f3890487
let fig = Figure(),
	dates = df_year[!, :dates],
    ax = Axis(fig[1,1];
        title="Time Series",
        xlabel="day",
        ylabel="Temperature (°C)")
	
    for ts in eachcol(select(df_year, Not(:dates)))[begin:32]
        lines!(ax, dates, ts)
    end

    fig
end

# ╔═╡ 9c69ee04-2dde-4c16-845c-4da70bfd4085
const M_ts = convert(Array{Float64}, reduce(hcat, eachcol(select(df_year, Not(:dates)))))

# ╔═╡ f6155294-9c55-4c99-b19d-2241c5291708
const n_samples = 1_000

# ╔═╡ ef36e6ad-6388-4547-9985-4a0c1b604283
const n_series = 100

# ╔═╡ 267d0c0d-c364-47d6-a652-0d5ec49127ba
const n_bins = 256

# ╔═╡ 61fc4161-d0d1-40dd-b96f-c65af38e5e63
@inline generate_columns(M::AbstractMatrix{<:Real}, n_series::Integer, n_samples::Integer) =
    let n_cols = size(M, 2)
    	map(_ -> M[:, randperm(n_cols)[begin:n_series]], 1:n_samples)
    end

# ╔═╡ 7671c085-97e8-440e-9c98-89e527fb4f8a
const cov_eigvals = mapreduce(TimeSeries.covariance_matrix_eigvals, hcat, generate_columns(M_ts, n_series, n_samples))

# ╔═╡ 71edb0f9-0f18-448b-bad7-4b68a73805dd
let fig = Figure(),
	λ = cov_eigvals |> vec |> filter(>(0)) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance eigenvalues Zipf plot",
        xlabel=L"i",
        ylabel=L"\lambda_i", yscale=log10)

    scatter!(ax, λ)

    fig
end

# ╔═╡ 77aa6109-6e43-4f61-b0c5-5d1da4d81e20
let fig = Figure(),
	λ = cov_eigvals |> filter(>=(0)) |> sort,
    ax = Axis(fig[1,1];
        title="Covariance matrix spectrum CDF",
        xlabel=L"\lambda_i", xscale=log10,
        ylabel=L"F(\lambda)",
		limits=(extrema(λ), (0, nothing)))

    scatterlines!(ax, λ, range(0, 1, length(λ)))

    fig
end

# ╔═╡ 153ce6bd-51a8-4539-9aa0-391f032e9f1a
let fig = Figure(),
	λ = cov_eigvals |> filter(>(1e-6)) |> sort,
	ax = Axis(fig[1,1], yscale=log10)

    hist!(ax, λ, bins=128, normalization=:pdf)

	fig
end

# ╔═╡ bc9b7623-f979-4eff-ad5d-e5c914e30c13
let fig = Figure(),
	λ = cov_eigvals |> filter(>(3e-6)) |> sort,
	hist = HistStats.LogHistogram(λ, n_bins),
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

# ╔═╡ b3459a5f-40e2-4536-8939-3ced6edb827e
let fig = Figure(),
	λ = filter(>=(0), cov_eigvals) |> sort,
	(mp, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(365, n_series)
	ax = Axis(fig[1,1];
        title="Covariance matrix spectrum",
        xlabel=L"\lambda_i",
        ylabel=L"\rho(\lambda)", yscale=Makie.log10)

    #barplot!(ax, x, hist.freqs, gap=0)
	hist!(ax, λ, bins=128, normalization=:pdf)

	x = 0:0.001:1000
	lines!(ax, x, mp.(x), color=:red)

	fig
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8f814715-3830-42be-85d4-1f29c2b0bffb
# ╟─b3dc1a6f-769e-43b3-ac16-f2d23e8cb71f
# ╟─8002063f-c2a7-40d8-b8c8-8c244ecec2e9
# ╟─3f0226c8-bce6-4c30-992c-6a2c61cd5dd1
# ╟─70c494fb-b3f8-432d-8948-ff1e690aeaa2
# ╟─9b4958f0-a21e-47a5-a246-9e05fe952f90
# ╟─b8212742-d884-4819-8f32-4fe9f3890487
# ╟─9c69ee04-2dde-4c16-845c-4da70bfd4085
# ╟─f6155294-9c55-4c99-b19d-2241c5291708
# ╟─ef36e6ad-6388-4547-9985-4a0c1b604283
# ╟─267d0c0d-c364-47d6-a652-0d5ec49127ba
# ╠═61fc4161-d0d1-40dd-b96f-c65af38e5e63
# ╠═7671c085-97e8-440e-9c98-89e527fb4f8a
# ╠═71edb0f9-0f18-448b-bad7-4b68a73805dd
# ╠═77aa6109-6e43-4f61-b0c5-5d1da4d81e20
# ╠═153ce6bd-51a8-4539-9aa0-391f032e9f1a
# ╠═bc9b7623-f979-4eff-ad5d-e5c914e30c13
# ╠═b3459a5f-40e2-4536-8939-3ced6edb827e
