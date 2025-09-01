module SeqFold

import TOML: parsefile
using Base: VersionNumber

const VERSION = let
    proj_path = joinpath(dirname(@__DIR__), "Project.toml")
    isfile(proj_path) ? VersionNumber(parsefile(proj_path)["version"]) : v"0.0.0"
end

"""
    SeqFold.jl, v$VERSION

Nucleic acid folding library

# Exports:
[`tm`](@ref), [`tm_cache`](@ref), [`gc_cache`](@ref), [`MeltingConditions`](@ref),
[`fold`](@ref), [`dg`](@ref), [`dg_cache`](@ref), [`dot_bracket`](@ref)

# Documentation:
https://phlaster.github.io/SeqFold.jl

$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "# Local path:\n" * pathof(SeqFold)) : "") 
"""
SeqFold

export tm, tm_cache, gc_cache, MeltingConditions
export fold, dg, dg_cache, dot_bracket
public complement

using Printf

include("utils.jl")
include("types.jl")
include("rna.jl")
include("dna.jl")
include("tm.jl")
include("fold.jl")

end # module