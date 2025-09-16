module SeqFold

"""
Nucleic acid folding and melting temperature calculation

Exports: [`tm`](@ref), [`MeltingConditions`](@ref), [`fold`](@ref), [`dg`](@ref), [`dot_bracket`](@ref)

Documentation: https://phlaster.github.io/SeqFold.jl/

$(isnothing(get(ENV, "CI", nothing)) ? ("\nLocal path: $(pathof(SeqFold))") : "") 
"""
SeqFold

export tm, tm_deg, MeltingConditions
export dg, fold, dot_bracket
public complement, Structure, tm_cache, gc_cache, dg_cache

using Printf

include("utils.jl")
include("types.jl")
include("rna.jl")
include("dna.jl")
include("tm.jl")
include("fold.jl")

end # module