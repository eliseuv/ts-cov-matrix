module DataIO

export
    # DataFile type
    DataFile,
    # Filename
    filename,
    # Parsing
    parse_extension, parse_filename,
    # Check parameters
    check_params,
    # Find datafiles
    list_datafiles,
    find_datafiles,
    # Datafiles vector
    filter_params,
    sort_datafiles,
    #= # Pickle =#
    #= load_pickle, =#
    #= load_pickle_gz, =#
    #= write_pickle_gz, =#
    # Rust ndarray
    load_ndarray,
    # Loading
    load

using GZip, JLD2, CBOR

include("Macros.jl")
using .Macros

# Default file extension
const DEFAULT_DATAFILE_EXT = ".cbor.gz"

#########################
# Parameters dictionary #
#########################

@doc raw""""
    merge_params(params::Union{Dict{String},Pair{String}}...)::Dict{String,Any}

Merge parameters pairs and dictionaries of parameters into a single dictionary.
"""
@inline merge_params(params::Union{Dict{String},Pair{String}}...)::Dict{String,Any} = reduce(merge, map(Dict, params); init=Dict{String,Any}())

@inline var_str((name, value)::Pair{String})::String =
    name * "=" * string(value)

@doc raw"""
    params_str(params...; sep="_")

Generate a string containing the 'name=value' of the parameters specified by the arguments `params`
(name-value pairs or dictionaries) sorted alphabetically and separated by the string `sep`.

# Example:
    ```julia
    julia> params_str(Dict("foo" => 0.5), "bar" => "baz")
    "bar=baz_foo=0.5"
    julia> params_str("foo" => 0.5, "bar" => "baz", sep=";")
    "bar=baz;foo=0.5"
    ```

"""
@inline params_str(params::AbstractVector{Pair{String,Any}}; sep::AbstractString="_")::String =
    join(map(var_str, sort(params, by=first)), sep)

@inline params_str(params::Dict{String,Any}; sep::AbstractString="_")::String =
    params_str([params...], sep=sep)

@inline params_str(params::Union{Dict{String},Pair{String}}...; sep::AbstractString="_")::String =
    params_str(merge_params(params...), sep=sep)

####################
# Filename parsing #
####################

@inline add_dot(ext::AbstractString) =
    if !startswith(ext, '.')
        '.' * ext
    else
        ext
    end

@doc raw"""
    filename(prefix::AbstractString, params...; sep::AbstractString="_", ext::AbstractString=DEFAULT_DATAFILE_EXT)

Generate a filename give an `prefix`, dictionaries of parameters, or pairs `params...` and a file extension `ext`.

Each parameter is written as `param_name=param_value` and separated by a `sep` string.

The dot `.` in the extension can be omitted: `ext=".csv"` and `ext="csv"` are equivalent.

The default file extension is `DEFAULT_DATAFILE_EXT`.
To create a file without extension, use either `ext=nothing` or `ext=""`.
"""
function filename(prefix::AbstractString, params...; sep::AbstractString="_", ext::AbstractString=DEFAULT_DATAFILE_EXT)
    # Prefix and parameters
    filename = prefix
    parameters = params_str(params..., sep=sep)
    if !isempty(parameters)
        filename *= sep * parameters
    end
    # Extension
    if !isnothing(ext) && ext != ""
        if ext[begin] == '.'
            filename *= ext
        else
            filename *= '.' * ext
        end
    end
    return filename
end

function parse_extension(path::AbstractString)
    (base, ext) = splitext(path)
    if ext == ".gz"
        (base, ext_second) = splitext(base)
        ext = ext_second * ext
    end
    return (base, ext)
end

@doc raw"""
    parse_filename(path::AbstractString; sep::AbstractString = "_")

Attempts to parse parameters in name of file given by `path` using `sep` as parameter separator.

It assumes the following pattern for the filename (using the default separator `"_"`):
    `SomePrefix_first_param=foo_second_param=42_third_param=3.14.ext`

Returns a tuple `(String, Dict{String,Any}, String)` containing:
    [1] Filename prefix
    [2] Dictionary with keys being the names of the parameters as symbols and the values the parsed parameter values
    [3] File extension

# Example:
    ```julia
    julia> parse_filename("/path/to/SomePrefix_first_param=foo_second_param=42_third_param=3.14.ext")
    ("SomePrefix", Dict("first_param" => "foo", "second_param" => 42, "third_param" => 3.14), ".ext")
    ```
"""
function parse_filename(path::AbstractString; sep::AbstractString="_")
    # Get filename and extension
    (filename, ext) = path |> basename |> parse_extension
    # Split name into chunks
    namechunks = split(filename, sep)
    # The first chunk is always the prefix
    prefix = popfirst!(namechunks)
    # Dictionary to store parsed parameter values
    params_dict = Dict{String,Any}()
    while length(namechunks) != 0
        param = popfirst!(namechunks)
        while !occursin("=", param) && length(namechunks) != 0
            param = param * sep * popfirst!(namechunks)
        end
        if occursin("=", param)
            (param_name, param_value) = split(param, "=")
        else
            break
        end
        # Try to infer type
        ParamType = infer_type_sized(param_value)
        if ParamType != Any
            # Type could be inferred, parse it
            params_dict[string(param_name)] = parse(ParamType, param_value)
        else
            # Type could not be inferred, keep it as String
            params_dict[string(param_name)] = param_value
        end
    end
    return (prefix, params_dict, ext)
end

##################
# DataPoint type #
##################

@doc raw"""
    DataPoint
"""
struct DataPoint
    prefix::AbstractString
    params::Dict{String,T} where {T}
end

#################
# DataFile type #
#################

@doc raw"""
    DataFile
"""
struct DataFile

    path::AbstractString
    dir::AbstractString
    datapoint::DataPoint
    ext::AbstractString

    # New datafile from parameters
    DataFile(dir::AbstractString, prefix::AbstractString, params...; ext::AbstractString=DEFAULT_DATAFILE_EXT, sep::AbstractString="_") =
        let path = joinpath(dir, filename(prefix, params...; sep=sep, ext=ext))
            datapoint = DataPoint(prefix, merge_params(params...))
            new(path, dir, datapoint, ext)
        end

    # New datafile from path
    DataFile(path::AbstractString; sep::AbstractString="_") =
        let path = abspath(path),
            (dir, filename) = splitdir(path),
            (prefix, params, ext) = parse_filename(filename; sep=sep)

            datapoint = DataPoint(prefix, params)

            new(path, dir, datapoint, ext)
        end
end

####################
# Check parameters #
####################

@doc raw"""
    check_params(params::Dict{String}, (key, value)::Pair{String,T}) where {T}

Checks if the prameters dictionary `params` has the key-value pair specified by `key => value`.
"""
@inline check_params(params::Dict{String}, (key, value)::Pair{String}) =
    haskey(params, key) && params[key] == value

@doc raw"""
    check_params(params::Dict{String}, (key, values)::Pair{String,<:AbstractVector{T}}) where {T}

Checks if the prameters dictionary `params` satisfies the any of the values specified by the pair `key => values`.
"""
@inline check_params(params::Dict{String}, (key, values)::Pair{String,<:AbstractVector}) =
    haskey(params, key) && any(map(value -> params[key] == value, values))

@doc raw"""
    check_params(params::Dict{String}, (key, predicate)::Pair{String,Function})

Checks if the prameters dictionary `params` satisfies the predicate specified by the pair `key => predicate`.
"""
@inline check_params(params::Dict{String}, (key, predicate)::Pair{String,<:Function}) =
    haskey(params, key) && predicate(params[key])

@doc raw"""
    check_params(params::Dict{String}, reqs::Pair{String}...)

Checks if the parameters dictionary `params` satisfies the values defined in the parameters pairs `reqs...`.
"""
@inline check_params(params::Dict{String}, reqs::Pair{String}...) =
    all(map(req -> check_params(params, req), reqs))

@doc raw"""
    check_params(params::Dict{String}, reqs::Dict{String})

Checks if the parameters dictionary `params` satisfies the values defined in the parameters requirements dictionary `reqs`.
"""
@inline check_params(params::Dict{String}, reqs::Dict{String}) =
    check_params(params, reqs...)

@doc raw"""
    check_params(params::Dict{String}, reqs...)

Checks if the parameters dictionary `params` satisfies the values defined in the parameters pairs or dictionaries `reqs...`.
"""
@inline check_params(params::Dict{String}, reqs...) = all(check_params(params, req) for req ∈ reqs)

@doc raw"""
    check_params(datafile::DataFile, reqs...)
"""
@inline check_params(datafile::DataFile, reqs...) = check_params(datafile.datapoint.params, reqs...)

##################
# Find datafiles #
##################

@inline list_datafiles(path::AbstractString; ext::AbstractString=DEFAULT_DATAFILE_EXT) =
    collect(Iterators.flatmap(walkdir(path)) do (root, _, filenames)
        map(filename -> joinpath(root, filename), filenames)
    end) |>
    filter(path -> parse_extension(path)[2] == add_dot(ext))

"""
    find_datafiles(path::AbstractString, prefix::AbstractString, reqs...; ext::AbstractString=DEFAULT_DATAFILE_EXT, sep::AbstractString="_")

Find data files in the directory `datadir` that have the satisfies the required parameters `reqs...`.
"""
@inline find_datafiles(path::AbstractString, prefix::AbstractString, reqs...; ext::AbstractString=DEFAULT_DATAFILE_EXT, sep::AbstractString="_") =
    collect(Iterators.flatmap(walkdir(path)) do (root, _, filenames)
        map(filename -> joinpath(root, filename), filenames)
    end) |>
    filter(path -> parse_extension(path)[2] == add_dot(ext)) |>
    paths -> map(path -> DataFile(path; sep=sep), paths) |>
             datafiles -> filter(datafile -> datafile.datapoint.prefix == prefix && check_params(datafile, reqs...), datafiles)

#####################
# Filter parameters #
#####################

@inline filter_params(datafiles::AbstractVector{DataFile}, reqs...) = filter(datafile -> check_params(datafile, reqs...), datafiles)

@inline filter_params(reqs...) = filter(datafile -> check_params(datafile, reqs...))

##################
# Sort datafiles #
##################

@inline get_params_values(datafiles::AbstractVector{DataFile}, param_name::String) = map(datafile -> datafile.datapoint.params[param_name], datafiles) |> sort |> unique
@inline function only_verbose(vec::AbstractVector)
    if (local n = length(vec)) > 1
        println("Error on `only`: length = $n")
        @show vec
    end
    return only(vec)
end

@inline sort_datafiles(datafiles::AbstractVector{DataFile}, param_name::String) =
    map(get_params_values(datafiles, param_name)) do param_value
        param_value => filter_params(datafiles, param_name => param_value) |> only_verbose
    end

@inline sort_datafiles(datafiles::AbstractVector{DataFile}, param_names::String...) =
    let param_name = param_names[begin]
        map(get_params_values(datafiles, param_name)) do param_value
            param_value => sort_datafiles(filter_params(datafiles, param_name => param_value), param_names[2:end]...)
        end
    end

@inline sort_datafiles(param_names::String...) =
    @inline f(datafiles::AbstractVector{DataFile}) = sort_datafiles(datafiles, param_names...)

#= #################### =#
#= # Pickle datafiles # =#
#= #################### =#
#==#
#= @pyinclude(joinpath(PROJECT_DIR, "data_io.py")) =#
#==#
#= @inline load_pickle = py"load_pickle" =#
#= @inline load_pickle_gz = py"load_pickle_gz" =#
#= @inline write_pickle_gz = py"write_pickle_gz" =#

################
# Rust ndarray #
################

function load_ndarray(dict::Dict)::Union{Nothing,Array{<:Real}}
    # Check keys
    for key ∈ ["v", "dim", "data"]
        if !haskey(dict, key)
            @error "Missing key in dict: $key"
            return nothing
        end
    end
    # Check version
    if dict["v"] != 1
        @error "Unsupported version!"
        return nothing
    end

    dims = tuple(reverse(dict["dim"])...)
    return reshape(dict["data"], dims)

end

##################
# Load datafiles #
##################

@inline load(datafile::DataFile) =
#= if datafile.ext == ".pickle" =#
#=     load_pickle(datafile.path) =#
#= elseif datafile.ext == ".pickle.gz" =#
#=     load_pickle_gz(datafile.path) =#
    if datafile.ext == ".cbor"
        open(datafile.path) do file
            read(file) |> CBOR.decode
        end
    elseif datafile.ext == ".cbor.gz"
        GZip.open(datafile.path) do file
            read(file) |> CBOR.decode
        end
    elseif datafile.ext == ".jld2"
        load_object(datafile.path)
    else
        println("Unrecognized file extension: $(datafile.ext)")
        nothing
    end

end
