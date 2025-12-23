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
    using Statistics, LinearAlgebra, JLD2, DataFrames, CairoMakie, ColorSchemes, LsqFit, NonlinearSolve
    DataIO = ingredients("./DataIO.jl").DataIO
    TimeSeries = ingredients("./TimeSeries.jl").TimeSeries
    HistStats = ingredients("./HistStats.jl").HistStats
end

# ╔═╡ 8b76c759-d2e3-449d-a4d1-7ae1fd6f0271
set_theme!(Theme(fontsize=24))

# ╔═╡ c86881af-46be-4787-a99d-6385c85c2099
@. belehradek_model(x, p) = p[1] * (x - p[2])^p[3]

# ╔═╡ 1339743f-08ab-4601-ad4a-ba15e197007b
begin
	const n_steps = 300
	const n_runs = 1000
end

# ╔═╡ 6b518705-7913-479b-a2ea-0a7ebfc17a38
datafiles = DataIO.find_datafiles("./data/contact_process/1d/all_active/eigvals/",
					  "ContactProcessDiffusion1DEigvals",
					  ext=".jld2",
					  "L" => 128) |> DataIO.sort_datafiles("diffusion", "rate")

# ╔═╡ a09a5092-925c-4523-88b5-e6fe1c1e1416
df_dict = map(datafiles) do (diffusion, datafiles_rate)
	diffusion => map(datafiles_rate) do (rate, datafile)
		eigvals_matrix = DataIO.load(datafile)
		hist = HistStats.Histogram(vec(eigvals_matrix), 128)
		eigvals_hist_mean = mean(hist)
		eigval_max_mean = mean(eigvals_matrix[end,:])
		(rate=rate,
		 eigvals_hist_mean=eigvals_hist_mean,
		 eigvals_var=var(eigvals_matrix),
		 eigval_max_mean=eigval_max_mean,
		 eigval_max_var=varm(eigvals_matrix[end,:], eigval_max_mean))
	end |> DataFrame
end |> Dict

# ╔═╡ c2dbe90e-6fd1-4c57-8755-a60297e36bf8
let fig = Figure(size=(800, 600)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(df_dict), categorical = true),
	ax = Axis(fig[1,1];
			  yscale=log10,
			  limits=(nothing, nothing))

	inflection_points = []
	for ((diffusion, df), color) in zip(sort(collect(df_dict), by=first), cmap)
		
		x = df[!, :rate]
		y = df[!, :eigval_max_var]
		
		lines!(ax,
			   x, y,
			   label="$(diffusion)", color=color)
		
	end

	Colorbar(fig[begin,end+1], limits=(0, 1), colormap=cmap,)

	fig
end

# ╔═╡ ea1c79a3-1006-45a3-8e2f-50fa7f8cf65c
@. maximum_model(x, p) = p[2]*(x-p[1])^2 + p[3] + p[4]*(x-p[1])^4

# ╔═╡ 74854305-8253-4a51-be7d-36dce63dcb0c
let fig = Figure(),
	ax = Axis(fig[1,1]),
	df = df_dict[0],
	x = range(2, 4, length=100)

	scatter!(ax, df[!, :rate], log10.(df[!, :eigval_max_var]))
	p0 = [3.2, -5, 0.65, 3]
	lines!(ax, x, maximum_model(x, p0))

	fig
end

# ╔═╡ 6e7302f6-8ed6-4862-836b-2f6390ba3522
eigval_max_var_plot, maxium_points = let fig = Figure(size=(1000, 600)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(df_dict), categorical = true),
	ax = Axis(fig[1,1];
			  title="(b) Maximum eigenvalue variance",
			  limits=((2, 4), nothing),
			  yscale=log10,
			  xlabel=L"\alpha",
			  ylabel=L"\text{var}(\lambda_\text{max})", ylabelrotation=0),
	estimates = [3.22, 3.1, 3, 2.92, 2.85, 2.79, 2.73, 2.69, 2.65, 2.6, 2.57, 2.54, 2.51, 2.48, 2.46, 2.44, 2.42, 2.4, 2.38, 2.36, 2.34],
	model = maximum_model

	max_points = []
	for ((diffusion, df), color, estimate) in zip(sort(collect(df_dict), by=first), cmap, estimates)
				
		lines!(ax,
			   df[!, :rate], df[!, :eigval_max_var],
			   label="$(diffusion)", color=color)

		L = 0.25
		df_fit = df[df.rate .- L .<= estimate .<= df.rate .+ L, :]
		x_fit = df_fit[!, :rate]
		y_fit = log10.(df_fit[!, :eigval_max_var])
		#p0 = [estimate, -2, 2, -1.4]
		p0 = [estimate, -5, 0.65, 3]
		fit = curve_fit(model, x_fit, y_fit, p0)
		alpha_crit = fit.param[1]
		y_crit = 10 ^ model(alpha_crit, fit.param)
		println(fit.param)

		#println("$(diffusion): $(errors)")
		push!(max_points, (diffusion=diffusion, alpha_crit=alpha_crit, y_crit=y_crit))
		
		x = range(extrema(df_fit.rate)..., length=100)
		lines!(ax, x, 10 .^ model(x, fit.param), color=:red, linestyle=:dash, alpha=1)
				
	end

	for p in max_points
		scatter!(ax, p.alpha_crit, p.y_crit, color=:red, markersize=8)
	end
	
	Colorbar(fig[begin,end+1], limits=(0, 1), colormap=cmap,
			 label=L"\gamma", labelrotation=0)

	fig, DataFrame(max_points)
end

# ╔═╡ 76a329fa-2eb1-4731-a9c7-0a1ed16cf293
eigval_max_var_plot

# ╔═╡ 85b60ef4-2546-4f3f-ae76-7f426f92d872
save("cp1d_eigval_max_var.pdf", eigval_max_var_plot)

# ╔═╡ 4d7b62fe-4a3a-420d-920d-3dbf649e150c
maxium_points

# ╔═╡ 9316ef0a-4bf3-4574-a4f2-0ed55e780221
let fig = Figure(size=(800, 600)),
	ax = Axis(fig[1,1];
			  title="Contact Process 1D with Diffusion",
			  xlabel=L"\gamma",
			  ylabel=L"\alpha_{c}", ylabelrotation=0,
			  limits=((0, 1), nothing))

	x = maxium_points[!, :diffusion]
	y = maxium_points[!, :alpha_crit]
	scatter!(ax, x, y, label="Maxima")

	fit = curve_fit(belehradek_model, x, y, [2.33,-0.16,-0.19])
	params = fit.param
	println(params)
	errors = standard_errors(fit)
	println(errors)
	
	x = range(extrema(x)..., length=100)
	lines!(ax, x, belehradek_model(x, params), label="Belehradek fit")

	axislegend(ax)
	
	fig

end

# ╔═╡ f6c719f5-b379-4242-a944-e3421466da2b
@. logistic_model(x, p) = p[4] + p[3] / (1 + exp(-p[2] * (x - p[1])))

# ╔═╡ 04ee620d-05c7-4358-81c5-95edbb1ba032
@. inflection_model(x, p) = p[3]*(x-p[1]) + p[4]*(x-p[1])^3 + p[2]

# ╔═╡ 9d25576f-d051-4aea-adc0-0c2edd5a05bd
let fig = Figure(),
	ax = Axis(fig[1,1]),
	df = df_dict[0],
	x = range(2, 4, length=100)

	scatter!(ax, df[!, :rate], df[!, :eigval_max_mean])
	p0 =[3.2, 50, -100, 100]
	lines!(ax, x, inflection_model(x, p0))

	fig
end

# ╔═╡ 846c3e22-a8c5-4ddc-be0d-d6bfd5afc59c
eigval_max_mean_plot, inflection_points = let fig = Figure(size=(1000, 600)),
	cbarPal = :viridis,
	cmap = cgrad(colorschemes[cbarPal], length(df_dict), categorical = true),
	ax = Axis(fig[1,1];
			  title="(a) Maximum eigenvalue average",
			  limits=((2, 4), nothing),
			  xlabel=L"\alpha",
			  ylabel=L"\langle\lambda_\text{max}\rangle", ylabelrotation=0),
	estimates = [3.22, 3.1, 3, 2.92, 2.85, 2.79, 2.73, 2.69, 2.65, 2.6, 2.57, 2.54, 2.51, 2.48, 2.46, 2.44, 2.42, 2.4, 2.38, 2.36, 2.34],
	model = logistic_model

	max_points = []
	for ((diffusion, df), color, estimate) in zip(sort(collect(df_dict), by=first), cmap, estimates)
				
		lines!(ax,
			   df[!, :rate], df[!, :eigval_max_mean],
			   label="$(diffusion)", color=color)

		L = 0.256
		df_fit = df[df.rate .- L .<= estimate .<= df.rate .+ L, :]
		x_fit = df_fit[!, :rate]
		y_fit = df_fit[!, :eigval_max_mean]
		p0 = [estimate, -5, 85, 12]
		#p0 =[estimate, 50, -100, 100]

		fit = curve_fit(model, x_fit, y_fit, p0)
		alpha_crit = fit.param[1]
		y_crit = model(alpha_crit, fit.param)
		#println(fit.param[1])

		#println("$(diffusion): $(errors)")
		push!(max_points, (diffusion=diffusion, alpha_crit=alpha_crit, y_crit=y_crit))
		
		x = range(extrema(df_fit.rate)..., length=100)
		lines!(ax, x, model(x, fit.param), color=:red, linestyle=:dash, alpha=1)
		
	end

	for p in max_points
		scatter!(ax, p.alpha_crit, p.y_crit, color=:red, markersize=8)
	end
	
	Colorbar(fig[begin,end+1], limits=(0, 1), colormap=cmap,
			 label=L"\gamma", labelrotation=0)

	fig, DataFrame(max_points)
end

# ╔═╡ 4e89f6ba-822c-49af-82b9-c954256c3f8e
eigval_max_mean_plot

# ╔═╡ 46b792e7-cefa-4068-a3f2-7decec4686a6
save("cp1d_eigval_max_mean.pdf", eigval_max_mean_plot)

# ╔═╡ 6acdaa7c-bf6a-4be0-860b-c8c8d18177d6
inflection_points

# ╔═╡ 89eb8d2f-5de3-4f18-9052-500f42b8ea96
let fig = Figure(size=(800, 600)),
	ax = Axis(fig[1,1];
			  title="Contact Process 1D with Diffusion",
			  xlabel=L"\gamma",
			  ylabel=L"\alpha_{c}", ylabelrotation=0,
			  limits=((0, 1), nothing))

	x = inflection_points[!, :diffusion]
	y = inflection_points[!, :alpha_crit]
	scatter!(ax, x, y, label="Inflection points")

	fit = curve_fit(belehradek_model, x, y, [2.33,-0.16,-0.19])
	params = fit.param
	println(params)
	errors = standard_errors(fit)
	println(errors)
	
	x = range(extrema(x)..., length=100)
	lines!(ax, x, belehradek_model(x, params), label="Belehradek fit")

	axislegend(ax)
	
	fig

end

# ╔═╡ Cell order:
# ╠═545aaffc-f9e6-11ef-3071-a1ea565f5cb0
# ╠═06f4e22b-9379-46a8-9477-824aa10f7bee
# ╠═8b76c759-d2e3-449d-a4d1-7ae1fd6f0271
# ╠═c86881af-46be-4787-a99d-6385c85c2099
# ╠═1339743f-08ab-4601-ad4a-ba15e197007b
# ╠═6b518705-7913-479b-a2ea-0a7ebfc17a38
# ╠═a09a5092-925c-4523-88b5-e6fe1c1e1416
# ╠═c2dbe90e-6fd1-4c57-8755-a60297e36bf8
# ╠═ea1c79a3-1006-45a3-8e2f-50fa7f8cf65c
# ╠═74854305-8253-4a51-be7d-36dce63dcb0c
# ╠═6e7302f6-8ed6-4862-836b-2f6390ba3522
# ╠═76a329fa-2eb1-4731-a9c7-0a1ed16cf293
# ╠═85b60ef4-2546-4f3f-ae76-7f426f92d872
# ╠═4d7b62fe-4a3a-420d-920d-3dbf649e150c
# ╠═9316ef0a-4bf3-4574-a4f2-0ed55e780221
# ╠═f6c719f5-b379-4242-a944-e3421466da2b
# ╠═04ee620d-05c7-4358-81c5-95edbb1ba032
# ╠═9d25576f-d051-4aea-adc0-0c2edd5a05bd
# ╠═846c3e22-a8c5-4ddc-be0d-d6bfd5afc59c
# ╠═4e89f6ba-822c-49af-82b9-c954256c3f8e
# ╠═46b792e7-cefa-4068-a3f2-7decec4686a6
# ╠═6acdaa7c-bf6a-4be0-860b-c8c8d18177d6
# ╠═89eb8d2f-5de3-4f18-9052-500f42b8ea96
