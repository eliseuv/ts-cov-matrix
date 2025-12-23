using DataFrames, CSV

include("TimeSeries.jl")
using .TimeSeries

const mc_steps = 365
const n_samples = 100
@show mc_steps, n_samples


(mp, (λ₋, λ₊)) = marchenko_pastur(mc_steps, n_samples)

@show λ₋, λ₊

const n_points = 256
lambda = range(λ₋, λ₊, n_points)
dist = map(mp, lambda)

df = DataFrame(lambda=lambda, dist=dist)
show(df)

output_datafile = "marchenko-pastur_steps=$(mc_steps)_samples=$(n_samples).csv"
@show output_datafile

CSV.write(output_datafile, df)

