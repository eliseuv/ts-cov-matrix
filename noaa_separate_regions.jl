using JLD2, DataFrames

df_ts = load_object("./data/noaa/1e5_uniform_points/fft_filter/L=0/time_series_L=0.jld2")
df_coords = load_object("./data/noaa/1e5_uniform_points/coords.jld2")

# Add index column
df_coords[!, :index] = 1:nrow(df_coords)

const latitude_cutoff = 40
@show latitude_cutoff

const output_dir = "./data/noaa/1e5_uniform_points/fft_filter/L=0/latitudes/"

# Low latitudes
df_ts_low = let indices = filter(:latitude => lat -> abs(lat) <= latitude_cutoff, df_coords)[!, :index]
    df_ts[!, ["dates", map(idx -> "x" * string(idx), indices)...]]
end
n_low = ncol(df_ts_low) - 1
println("Low latitudes points: $(n_low)")
save_object(joinpath(output_dir, "lat_low_time_series.jld2"), df_ts_low)

# High latitudes
df_ts_high = let indices = filter(:latitude => lat -> abs(lat) > latitude_cutoff, df_coords)[!, :index]
    df_ts[!, ["dates", map(idx -> "x" * string(idx), indices)...]]
end
n_high = ncol(df_ts_high) - 1
println("High latitudes points: $(n_high)")
save_object(joinpath(output_dir, "lat_high_time_series.jld2"), df_ts_high)

println("Total points: $(n_low + n_high)")
