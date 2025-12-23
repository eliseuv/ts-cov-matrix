# Calculate eigenvalues of the associated covariance matrix for a given time series matrix.

using LinearAlgebra, JLD2

include("Utils.jl")
using .Utils
include("DataIO.jl")
using .DataIO
include("TimeSeries.jl")
using .TimeSeries

@inline get_covariance_matrix_eigvals(matrix::AbstractMatrix{<:Real}, n_runs::Integer) =
    let (n_steps, n_samples_total) = size(matrix)
        @assert n_samples_total % n_runs == 0 "The number of runs must divide the number of samples"
        n_samples = n_samples_total ÷ n_runs
        dropdims(mapslices(eigvals ∘ covariance_matrix, reshape(matrix, (n_steps, n_samples, n_runs)), dims=[1, 2]), dims=2)
    end

@inline get_intermediate_path(path::AbstractString, dir::AbstractString)::String = joinpath(dirname(path), dir, basename(path))

@inline get_output_datafile(input::DataFile, n_steps::Integer, n_runs::Integer)::DataFile =
    let dir = joinpath(input.dir, "eigvals"),
        n_samples_total = input.datapoint.params["n_samples"],
        prefix = input.datapoint.prefix * "Eigvals"

        @assert n_samples_total % n_runs == 0 "The number of runs must divide the number of samples"
        n_samples = n_samples_total ÷ n_runs

        DataFile(dir, prefix, merge(input.datapoint.params, Dict("n_steps" => n_steps, "n_samples" => n_samples, "n_runs" => n_runs)), ext=".jld2")
    end

const n_runs = 1000
const n_steps = 500

# Datafiles
const datafiles = let paths = map(abspath, ARGS),
    input_datafiles = map(DataFile, paths) |> filter_params("rate" => x -> 2 <= x <= 4),
    output_datafiles = map(idf -> get_output_datafile(idf, n_steps, n_runs), input_datafiles)

    filter(x -> !isfile(x[2].path), collect(zip(input_datafiles, output_datafiles)))
end

const n = length(datafiles)
const field = "time_series_matrix"

@info "Processing $n datafiles"

for (i, (input_datafile, output_datafile)) ∈ enumerate(datafiles)

    @info "$i/$n"
    @show input_datafile.path
    @show output_datafile.path
    mkpath(dirname(output_datafile.path))

    @info "Loading datafile..."
    @info input_datafile.path
    data = try
        DataIO.load(input_datafile)
    catch e
        @error "Unable to load datafile\n$(e)"
        mv(input_datafile.path, get_intermediate_path(input_datafile.path, "error"))
        continue
    end
    repl_print(data)
    L = data["args"]["length"]
    m_ts = load_ndarray(data[field])[begin:n_steps+1, :] ./ L
    repl_print(m_ts)
    @info "Calculating..."
    normalize_ts!(m_ts)
    corr_mat_eigvals = get_covariance_matrix_eigvals(m_ts, n_runs)
    repl_print(corr_mat_eigvals)

    @info "Writing..."
    @info output_datafile.path
    save_object(output_datafile.path, corr_mat_eigvals)

end
