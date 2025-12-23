using JLD2, DataFrames, Dates, LsqFit

@doc """
 	approx_index(x::AbstractVector{<:Real}, x₀::Real)

 Find index `i` for which `x[i] ≈ x₀`.

 """
@inline approx_index(x::AbstractVector{<:Real}, x₀::Real) =
    if x[begin] > x₀
        findfirst(<=(x₀), x)
    else
        findfirst(>=(x₀), x)
    end

@doc """
 	init_params(x::AbstractVector{<:Real})

 Estimate initial parameters for the model fit.
 """
function init_params(x::AbstractVector{<:Real})
    # Offset
    (xₘᵢₙ, xₘₐₓ) = extrema(x)
    B = (xₘᵢₙ + xₘₐₓ) / 2

    # Amplitude
    A = (xₘₐₓ - xₘᵢₙ) / 2

    # Phase
    φ = approx_index(x, B) * 2π / 365

    return [B, A, φ]
end

model(t, p) = p[1] .+ p[2] * sin.(((2π / 365) * t) .+ p[3])

function remove_sin(x::AbstractVector{<:Real})
    t = 1:length(x)
    u₀ = init_params(x)
    fit = curve_fit(model, t, x, u₀)

    return x .- model(t, fit.param)
end

@info "Loading..."
df = load_object("./data/noaa/1e5_uniform_points/sst_time_series.jld2")

@info "Fitting..."
select!(df, :dates, Not(:dates) .=> remove_sin .=> Not(:dates))

@info "Saving..."
save_object("./data/noaa/1e5_uniform_points/sst_time_series_sin_removed.jld2", df)
