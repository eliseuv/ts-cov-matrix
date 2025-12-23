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
    using PlutoUI, Statistics, JLD2, Dates, DataFrames, FFTW, Plots, ColorSchemes, LaTeXStrings
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8da30af3-67dc-44f0-958e-57fdf85fa5e9
const df_ts = filter(:dates => d -> year(d) ∈ 1982:2024, load_object("./data/noaa/1e5_uniform_points/time_series.jld2"))

# ╔═╡ 96006342-7add-44eb-907e-eb5199fadb2b
const name = "x6"

# ╔═╡ 171f7d57-8c47-4d44-8b55-88cac95ba8ca
const x = df_ts[!, name]

# ╔═╡ 5a989c48-397d-472d-80cf-75ce4695eba9
# Number of samples
const Nₛ = length(x)

# ╔═╡ 5be2e3c5-7f62-4150-b741-0de428db1ec6
# Sampling rate
# Daily => Average number of days per year
const fₛ = 365.25

# ╔═╡ ae3c5797-9b9a-4ac8-96dd-57d3e9d42978
nrow(df_ts)/fₛ

# ╔═╡ 4b119a5c-7829-4fca-89a7-59ba1f97c335
# Time step
const Tₛ = 1 / fₛ

# ╔═╡ 7bf2c02e-4b76-4816-8392-69a1a898ddc7
# Time coordinates
const t =
	let t₀ = 0
		range(t₀, length=Nₛ, step=Tₛ)
	end

# ╔═╡ 4627933f-48a0-4697-a978-ad19328901aa
plot(t, x)

# ╔═╡ 05fce7ad-bb19-4a1c-9cb3-c98ac001d0c9
# Frequencies
const ω = rfftfreq(Nₛ, fₛ)

# ╔═╡ fda93894-ad6a-4b5f-a990-c0b08ada45f6
F = rfft(x)

# ╔═╡ 75937267-9177-48b8-a14f-5780e6531d0a
@bind ωₘₐₓ Slider(1:44; default=4)

# ╔═╡ 87aa3f06-9e54-4cce-968d-ffd6a121e442
let F = rfft(x),
	φ = map(angle, F),
	Fa = abs.(F),
	cmap = :twilight

	plt = plot(ω, Fa,
		 linewidth=2,
		 yscale=:log10,
		 xlim=(0, ωₘₐₓ),
		 label=nothing,
		 title="Time Series DFT")
	
	scatter!(ω, Fa,
			 markersize=3,
			 zcolor=φ, colorbar=true, c=cmap,
			 colorbar_title="φ", colorbar_titlefontrotation=π,
			 label=nothing)

	savefig("./plots/noaa/fft/fft_$(name).pdf")

	plt
	
end

# ╔═╡ 92ec432b-3ad9-42b3-88fc-9d10056ef6ed
@inline approx_idx(ran, val) =
	round(Int, (val - first(ran))/step(ran) + 1)

# ╔═╡ da86b125-c059-47f8-910d-ca600e5798b5
@inline band_stop_filter!(F, ω, ω₀, L=0) =
	let i = approx_idx(ω, ω₀)

		F[i-L:i+L] .= 0
		
	end

# ╔═╡ db93bc60-3dd3-4f2b-a85b-07f824417926
@inline function apply_filter!(F, ω)

	nₘₐₓ = floor(Int, ω[end])
	
	for n ∈ 1:nₘₐₓ
		band_stop_filter!(F, ω, n, 0)
	end
	
end

# ╔═╡ fb23790b-ec53-4168-b556-d387134ca037
let F = rfft(x),
	cmap = :twilight

	apply_filter!(F, ω)
	φ = map(angle, F)
		
	plot(ω, abs.(F),
		 linewidth=2,
		 yscale=:log10,
		 xlim=(0, ωₘₐₓ), ylim=(1e1, 3e4),
		 label=nothing,
		 title="Time Series DFT Filtered")
	
	scatter!(ω, abs.(F),
			 markersize=3,
			 zcolor=φ, colorbar=true, c=cmap,
			 colorbar_title="φ", colorbar_titlefontrotation=π,
			 label=nothing)

end

# ╔═╡ 5a8b78d9-c2d8-4159-b3a7-d17f78be2989
let F = rfft(x),
	ylim = (1e1, 1e5)
			
	original = plot(ω, abs.(F),
		 linewidth=2,
		 yscale=:log10, ylim=ylim,
		 label=nothing,
		 title="Original")
	
	scatter!(ω, abs.(F),
			 markersize=3,
			 label=nothing)

	apply_filter!(F, ω)
				
	filtered = plot(ω, abs.(F),
		 linewidth=2,
		 xlabel=L"\omega",
		 yscale=:log10, ylim=ylim,
		 label=nothing,
		 title="Filtered")
	
	scatter!(ω, abs.(F),
			 markersize=3,
			 label=nothing)

	plt = plot(original, filtered,
		 layout=(2, 1),
		 link=:x, xlim=(0, ωₘₐₓ))
	
	savefig(plt, "./plots/noaa/fft/fft_$(name).pdf")

	plt

end

# ╔═╡ ef491f45-1740-4f22-932e-2cbba59eda15
let F = rfft(x)

	apply_filter!(F, ω)

	x′ = irfft(F, Nₛ)

	orignal = plot(t, x, title="Original", label=nothing)
	removed = plot(t, x .- x′, title="Seasonality", label=nothing,
				   ylabel="Temperature")
	filtered = plot(t, x′, title="Result", label=nothing,
				    xlabel="Year")

	plt = plot(orignal, removed, filtered,
		 layout=(3, 1),
		 link=:x, xlim=extrema(t))

	savefig(plt, "./plots/noaa/fft/season_$(name).pdf")

	plt
end

# ╔═╡ 07e3d9d5-fd23-4742-8508-c11a7eb14bd3
@inline high_pass(s) =
	let τ = 2,
		k = 1/(2π*τ)
	(k*s) / (1 + k*s)
	end

# ╔═╡ 45ad53b5-196f-4365-a89c-b18c2a7de26b
@inline harmonic_attenuation!(F, ω) =
	let nₘₐₓ = floor(Int, ω[end])
	
	for n ∈ 1:nₘₐₓ
		i = approx_idx(ω, n)
		F[i] = high_pass(n) * F[i]
	end
		
end

# ╔═╡ 1516062e-6c24-4f14-82da-6f04386f1e9a
let F = rfft(x)

	orignal = plot(ω, abs.(F),
				   title="Original", label=nothing,
				   yscale=:log10)
	scatter!(ω, abs.(F), markersize=1.5, label=nothing)

	filter = plot(ω, high_pass,
				  title="High-pass filter", label=nothing,
				  yscale=:log10, ylims=(1e-4, 1))

	
	harmonic_attenuation!(F, ω)

	result =  plot(ω, abs.(F),
				   title="Result", label=nothing,
				   yscale=:log10)
	scatter!(ω, abs.(F), markersize=1.5, label=nothing)

	plot(orignal, filter, result,
		 layout=(3, 1),
		 link=:x, xlims=(0, ωₘₐₓ))
end

# ╔═╡ 5f55745a-439e-4944-879a-8f5d78bde121
let F = rfft(x),
	cmap = :twilight

	harmonic_attenuation!(F, ω)
	φ = map(angle, F)

	plot(ω, abs.(F),
		 linewidth=2,
		 yscale=:log10,
		 xlim=(0, ωₘₐₓ),
		 label=nothing)
	
	scatter!(ω, abs.(F),
			 markersize=3,
			 zcolor=φ, colorbar=true, c=cmap,
			 colorbar_title="φ", colorbar_titlefontrotation=π,
			 label=nothing)

end

# ╔═╡ bccea984-58b1-401e-85d3-ec88c98d0c6d
let F = rfft(x)

	harmonic_attenuation!(F, ω)

	x′ = irfft(F, Nₛ)

	orignal = plot(t, x, title="Original", label=nothing)
	removed = plot(t, x .- x′, title="Seasonality", label=nothing,
				   ylabel="Temperature")
	filtered = plot(t, x′, title="Result", label=nothing,
				    xlabel="Year")

	plot(orignal, removed, filtered,
		 layout=(3, 1),
		 link=:x, xlim=extrema(t))
end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8da30af3-67dc-44f0-958e-57fdf85fa5e9
# ╠═96006342-7add-44eb-907e-eb5199fadb2b
# ╠═171f7d57-8c47-4d44-8b55-88cac95ba8ca
# ╟─5a989c48-397d-472d-80cf-75ce4695eba9
# ╠═5be2e3c5-7f62-4150-b741-0de428db1ec6
# ╠═ae3c5797-9b9a-4ac8-96dd-57d3e9d42978
# ╟─4b119a5c-7829-4fca-89a7-59ba1f97c335
# ╟─7bf2c02e-4b76-4816-8392-69a1a898ddc7
# ╠═4627933f-48a0-4697-a978-ad19328901aa
# ╟─05fce7ad-bb19-4a1c-9cb3-c98ac001d0c9
# ╟─fda93894-ad6a-4b5f-a990-c0b08ada45f6
# ╟─75937267-9177-48b8-a14f-5780e6531d0a
# ╟─87aa3f06-9e54-4cce-968d-ffd6a121e442
# ╟─92ec432b-3ad9-42b3-88fc-9d10056ef6ed
# ╟─da86b125-c059-47f8-910d-ca600e5798b5
# ╟─db93bc60-3dd3-4f2b-a85b-07f824417926
# ╠═fb23790b-ec53-4168-b556-d387134ca037
# ╠═5a8b78d9-c2d8-4159-b3a7-d17f78be2989
# ╠═ef491f45-1740-4f22-932e-2cbba59eda15
# ╟─07e3d9d5-fd23-4742-8508-c11a7eb14bd3
# ╟─45ad53b5-196f-4365-a89c-b18c2a7de26b
# ╟─1516062e-6c24-4f14-82da-6f04386f1e9a
# ╟─5f55745a-439e-4944-879a-8f5d78bde121
# ╟─bccea984-58b1-401e-85d3-ec88c98d0c6d
