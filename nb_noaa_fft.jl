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
    using Statistics, JLD2, Dates, DataFrames, Plots, FFTW
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8da30af3-67dc-44f0-958e-57fdf85fa5e9
const df_ts = load_object("./data/noaa/1e5_uniform_points/sst_time_series.jld2")

# ╔═╡ 171f7d57-8c47-4d44-8b55-88cac95ba8ca
x = df_ts[!, :x5]

# ╔═╡ 878730d2-e4c7-4760-bc83-3079f1ef6ce1
begin
	# Sample length
	const N = length(x)
	# Sample rate
	const fs = 365
	# Sample spacing
	const Ts = 1 / fs
	# Time coordinates
	const t =
		let t₀ = 0,
			tₘₐₓ = t₀ + (N-1) * Ts
				0:Ts:tₘₐₓ
		end
end

# ╔═╡ 3f8fa22d-04ea-49ec-9fcf-03093f226eac
time_idx(tᵥₐₗ::Real) = round(Int, tᵥₐₗ / step(t))

# ╔═╡ e0ad2adf-d296-470f-86cf-3f27f2b041a5
time_idx(2)

# ╔═╡ 576b87c5-e966-4030-a9bc-a85fd9683457
freq_idx(ω::Real)

# ╔═╡ 4627933f-48a0-4697-a978-ad19328901aa
plot(t, x)

# ╔═╡ ec145da1-f258-4a65-8274-e6c4c67400ae
freqs = rfftfreq(length(x), fs)

# ╔═╡ b8ece65c-adc8-4ac6-ac6e-474f60c8c9a6
freq_index(freqs::AbstractVector{<:Real}, ω::Real) =
	round(Int, (ω - freqs[begin]) / step(freqs)) + 1

# ╔═╡ 36857e91-02f7-4138-8a74-0dcc0ed26511


# ╔═╡ 52fe7dcd-a083-43ba-b520-8209879c0bdb
F = let F = rfft(x),
		i_const = freq_index(freqs, 0),
		i_year = freq_index(freqs, 1),
		d = 5

	#F[i_const] = 0
	#F[i_year-d:i_year+d] .= 0

	F
end

# ╔═╡ 675e65f7-06f4-4882-abc8-2013c6e7f623
plot(freqs, abs.(F), xlims=(0, 15), ylims=(8e0, 1e5), yscale=:log10)

# ╔═╡ 77109506-7d39-41f2-8b9d-7233fc97b3fe
x′ = irfft(F, N)

# ╔═╡ 8e803d50-4ec2-4e12-bc71-a31685d1a529
plot(t, [x .- mean(x), x′])

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8da30af3-67dc-44f0-958e-57fdf85fa5e9
# ╠═171f7d57-8c47-4d44-8b55-88cac95ba8ca
# ╠═878730d2-e4c7-4760-bc83-3079f1ef6ce1
# ╠═3f8fa22d-04ea-49ec-9fcf-03093f226eac
# ╠═e0ad2adf-d296-470f-86cf-3f27f2b041a5
# ╠═576b87c5-e966-4030-a9bc-a85fd9683457
# ╠═4627933f-48a0-4697-a978-ad19328901aa
# ╠═ec145da1-f258-4a65-8274-e6c4c67400ae
# ╠═b8ece65c-adc8-4ac6-ac6e-474f60c8c9a6
# ╠═36857e91-02f7-4138-8a74-0dcc0ed26511
# ╠═52fe7dcd-a083-43ba-b520-8209879c0bdb
# ╠═675e65f7-06f4-4882-abc8-2013c6e7f623
# ╠═77109506-7d39-41f2-8b9d-7233fc97b3fe
# ╠═8e803d50-4ec2-4e12-bc71-a31685d1a529
