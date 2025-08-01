using JLD2, DataFrames, Dates, Statistics

get_seasonality(x::AbstractVector{<:Real}, period::Integer) =
    let n = length(x)

        map(1:period) do i
            mean(x[i:period:n])
        end

    end

remove_seasonality(x::AbstractVector{<:Real}, period::Integer) =
    let n = length(x)

        x .- repeat(get_seasonality(x, period), n ÷ period + 1)[begin:n]

    end

@info "Loading..."
df = load_object("./data/noaa/1e5_uniform_points/sst_time_series.jld2")

@info "Fitting..."
for col ∈ names(df) |> filter(!=("dates"))
    df[!, col] = remove_seasonality(df[!, col], 365)
end

@info "Saving..."
save_object("./data/noaa/1e5_uniform_points/sst_time_series_seasonality_removed.jld2", df)
