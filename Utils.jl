@doc raw"""
    Auxiliary functions

"""
module Utils

export
    repl_print,
    nest_map, to_dict_nested

@doc raw"""
    repl_show(x...)

Print the entities `x...` to standard output while in a script in the same way it prints in a REPL session.
"""
function repl_print(x...)
    show(IOContext(stdout, :limit => true), "text/plain", x...)
    println()
end

@doc raw"""
    nest_map(f, arr)

Map over nested array.
"""
@inline nest_map(f::Function, arr::AbstractArray{T}) where {T} = map(f, arr)

@inline nest_map(f::Function, arr::AbstractArray{T}) where {T<:AbstractArray} =
    map(arr_inner -> nest_map(f, arr_inner), arr)

@inline nest_map(f::Function, arr::AbstractArray{Pair{K,V}}) where {K,V} =
    map(arr) do (key, value)
        key => f(value)
    end

@inline nest_map(f::Function, arr::AbstractArray{Pair{K,V}}) where {K,V<:AbstractArray} =
    map(arr) do (key, pairs_arr)
        key => nest_map(f, pairs_arr)
    end

@doc raw"""
    to_nest_dict(pairs::AbstractVector{T}) where {T<:Union{Pair,Pair{K,AbstractVector{T}}},K}

Converts nested `Pair` structure to nested `Dict` structure.

# Example:

    ```julia
    [foo => ['a' => 0, ...],
     bar => [42 => [:kek => "lol", ...], ...], ...]
    ```

    becomes

    ```julia
    Dict([foo => Dict(['a' => 0, ...]),
          bar => Dict([42 => Dict([:kek => "lol", ...]), ...]), ...])
    ```
"""
@inline to_nest_dict(pairs::Vector{Pair{K,V}}) where {K,V} =
    Dict(pairs)

@inline to_dict_nested(pairs::Vector{Pair{K,Vector{Pair{U,V}}}}) where {K,U,V} =
    map(pairs) do (key, inner_pairs)
        key => to_dict_nested(inner_pairs)
    end |> Dict

end
