@doc raw"""
    Metaprogramming and macros
"""
module Macros

export
    # Variables to dictionary
    @varsdict,
    # Quote strings
    @dq_str, @sq_str, @p_str,
    NUMBER_REGEXES,
    # Infer variable type from string
    infer_type, infer_type_sized,
    # Parse number from string
    parse

@doc raw"""
    @varsdict(vars::Symbol...)

Returns a dictionary with variable names as keys and its values.

# Example
    ```julia
    julia> foo=5; bar=2.3; baz="spam";
    julia> @varsdict foo bar baz
    Dict{String, Any} with 3 entries:
      "bar" => 2.3
      "baz" => "spam"
      "foo" => 5
    ```
"""
macro varsdict(vars::Symbol...)
    dict = Expr(:call, :Dict)
    for var in vars
        push!(dict.args, :($(string(var)) => $(esc(var))))
    end
    return dict
end

# Quote macros (Suggested by https://stackoverflow.com/a/20483464)

@doc raw"""
Double quotes quote string

# Usage:
    dq"with double quotes" -> "with double quotes"
"""
macro dq_str(s)
    "\"" * s * "\""
end

@doc raw"""
Simple quotes quote string

# Usage:
    sq"with simple quotes" -> 'with simple quotes'
"""
macro sq_str(s)
    "'" * s * "'"
end

@doc raw"""
    p"quote"

p"quote" macro

Since regular strings do not allow for escape chars and the r"quote" macro
compiles the regular expression when defined, this intermediate quote allows
us to store parts of regular expressions that can later be combined and compiled.

# Example:
    ```julia
    julia> pq1 = p"\\D+"
    julia> pq2 = p"\\d+"
    julia> re = Regex(pq1 * p"=" * pq2)
    julia> re
    r"\\D+=\\d+"
    ```
"""
macro p_str(s)
    s
end

# Regular expressions that match certain types
@doc raw"""
    NUMBER_REGEXES

Dictionary that associates each numerical type with the corresponding regular expression used to infer it.
"""
const NUMBER_REGEXES = Dict{Type,String}(
    # Real numbers
    Integer => p"[-+]?[0-9]*([eE][+]?[0-9]+)?",
    AbstractFloat => p"[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9]+)?",
    # Complex numbers
    Complex{Integer} => p"[-+]?[0-9]*([eE][+]?[0-9]+)?([ ]*)?[-+]([ ]*)?[0-9]*([eE][+]?[0-9]+)?([ ]*)[ij]",
    Complex{AbstractFloat} => p"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?([ ]*)?[-+]([ ]*)?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?([ ]*)[ij]",
)

@doc raw"""
    infer_type(value::AbstractString)::Type

Try to infer the type of the numerical value in the string `value`.

If no type could be inferred, returns `Any`.
"""
function infer_type(value_str::AbstractString)::Type
    # Test numerical types in a given order
    for Type in (Integer, AbstractFloat, Complex{Integer}, Complex{AbstractFloat})
        re = Regex(p"^" * NUMBER_REGEXES[Type] * p"$")
        if occursin(re, strip(value_str))
            return Type
        end
    end
    # If no numerical type could be inferred, return `Any`
    return Any
end

@inline infer_type_sized(value_str::AbstractString)::Type =
    let InferredType = infer_type(value_str)
        if InferredType == Integer
            Int64
        elseif InferredType == AbstractFloat
            Float64
        elseif InferredType == Complex{Integer}
            Complex{Int64}
        elseif InferredType == Complex{AbstractFloat}
            Complex{Float64}
        else
            InferredType
        end
    end

end
