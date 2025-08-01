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
    using Statistics, Dates, JLD2, DataFrames, StatsBase, CairoMakie
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 1f087673-fac5-4c28-ad0d-ceb1816bba36
const df_sst = load("./data/noaa/random_points/uniform.jld2", "df_sst")

# ╔═╡ d2f5e286-26b4-4816-9875-1645db4b647a
@inline persistence(x) =
    map(1:length(x)-1) do i
        if x[i+1] > x[i]
            findfirst(<=(x[i]), x[i+1:end])
        else
            findfirst(>=(x[i]), x[i+1:end])
        end
    end

# ╔═╡ 1123032b-b0d1-4ee8-a87a-e5f810206b35
const persistences = reduce(hcat, map(persistence, eachcol(df_sst[!,Not(:dates)])))

# ╔═╡ 2adebaee-5048-4fb6-ac8d-7010d184a391
DataFrame(persistences, :auto)

# ╔═╡ 64b87035-9e0b-4bd9-befe-94dea20f5dda
per_filtered = map(xs -> [x for x ∈ xs if !isnothing(x)], eachcol(persistences))

# ╔═╡ 4de58538-86db-4a89-8984-5f2152faf5bd
hist(persistences[:,1]; bins=128, axis=(;yscale=log10))

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═1f087673-fac5-4c28-ad0d-ceb1816bba36
# ╠═d2f5e286-26b4-4816-9875-1645db4b647a
# ╠═1123032b-b0d1-4ee8-a87a-e5f810206b35
# ╠═2adebaee-5048-4fb6-ac8d-7010d184a391
# ╠═64b87035-9e0b-4bd9-befe-94dea20f5dda
# ╠═4de58538-86db-4a89-8984-5f2152faf5bd
