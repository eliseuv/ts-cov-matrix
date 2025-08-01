### A Pluto.jl notebook ###
# v0.20.8

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
    using CairoMakie
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ d67b398a-d7ce-4a9a-ab17-7d119043c9d4
begin
	const n_steps = 365
	const n_series = 100
	(f, (λ₋, λ₊)) = TimeSeries.marchenko_pastur(n_steps, n_series)
end

# ╔═╡ be928c71-3586-4fe2-a532-4e4cf95c2906
let fig = Figure(),
	ax = Axis(fig[1,1];
			 title="Marchenko-Pastur",
			 xlabel=L"\lambda",),
	x = range(λ₋, λ₊, 512)

	lines!(ax, x, f.(x))

	fig
end

# ╔═╡ Cell order:
# ╟─545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═d67b398a-d7ce-4a9a-ab17-7d119043c9d4
# ╠═be928c71-3586-4fe2-a532-4e4cf95c2906
