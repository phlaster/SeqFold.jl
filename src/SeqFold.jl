module SeqFold

import TOML: parsefile
using Base: VersionNumber

const VERSION = let
    proj_path = joinpath(dirname(@__DIR__), "Project.toml")
    isfile(proj_path) ? VersionNumber(parsefile(proj_path)["version"]) : v"0.0.0"
end

"""
    SeqFold.jl, v$VERSION

Nucleic acid folding and thermodynamics library implementing the Zuker and Stiegler (1981) 
folding algorithm and nearest-neighbor thermodynamics for DNA/RNA analysis.

This Julia implementation is based on the Python library `seqfold` (https://github.com/Lattice-Automation/seqfold), 
with significant improvements:
- Fixed numerous bugs in the original `tm` implementation
- Added fine-grained control over ion concentrations for Tm calculations while preserving 
  two standard presets (`:pcr` and `:std`)
- Enhanced validation of buffer conditions to prevent physically impossible scenarios
- Verified against Biopython's `Bio.SeqUtils.MeltingTemp.Tm_NN` for melting temperature accuracy
- Maintained identical folding results compared to the original `seqfold` implementation

The folding algorithm strictly follows the Zuker approach with energy parameters from 
Breslauer et al. (1986) for DNA and Turner 2004 for RNA. Tm calculations incorporate 
salt corrections from Owczarzy et al. (2008) with support for Mg²⁺, Na⁺, K⁺, Tris, and dNTPs.

# Exports:
[`tm`](@ref), [`MeltingConditions`](@ref),
[`fold`](@ref), [`dg`](@ref), [`dot_bracket`](@ref)

$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "# Local path:\n" * pathof(SeqFold)) : "") 
"""
SeqFold

export tm, MeltingConditions
export dg, fold, dot_bracket
public complement, Structure, dg_cache, tm_cache, gc_cache

using Printf

include("utils.jl")
include("types.jl")
include("rna.jl")
include("dna.jl")
include("tm.jl")
include("fold.jl")

end # module