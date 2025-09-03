const Comp = Dict{Char, Char}
const Cache = Vector{Vector{Float64}}
const MultiBranch = NTuple{4, Float64}
const BpEnergy = Dict{String, Tuple{Float64, Float64}}
const LoopEnergy = NTuple{30, Tuple{Float64, Float64}}

mutable struct Energies
    BULGE_LOOPS::LoopEnergy
    COMPLEMENT::Comp
    DE::BpEnergy
    HAIRPIN_LOOPS::LoopEnergy
    MULTIBRANCH::MultiBranch
    INTERNAL_LOOPS::LoopEnergy
    INTERNAL_MM::BpEnergy
    NN::BpEnergy
    TERMINAL_MM::BpEnergy
    TRI_TETRA_LOOPS::Union{Nothing, BpEnergy}
    
    function Energies(
        bulge_loops::LoopEnergy,
        complement::Comp,
        de::BpEnergy,
        hairpin_loops::LoopEnergy,
        multibranch::MultiBranch,
        internal_loops::LoopEnergy,
        internal_mm::BpEnergy,
        nn::BpEnergy,
        terminal_mm::BpEnergy,
        tri_tetra_loops::Union{Nothing, BpEnergy}=nothing
    )
        return new(
            bulge_loops,
            complement,
            de,
            hairpin_loops,
            multibranch,
            internal_loops,
            internal_mm,
            nn,
            terminal_mm,
            tri_tetra_loops
        )
    end
end

"""
    Structure(e, desc, ij)

Represents a structural element within a predicted nucleic acid secondary structure.

A `Structure` object encapsulates the energetic contribution (`e`), a descriptive
string (`desc`) detailing the type of structural element (e.g., hairpin, stack,
bulge, multibranch), and the base-pairing information (`ij`).

This type is primarily used internally by the folding algorithm (see [`fold`](@ref))
to represent and cache the energetics of different structural motifs. It can also
be used to inspect the detailed components of a folded structure.

# Fields

- `e::Float64`: The free energy contribution (in kcal/mol) of this structural element.
  A value of `-Inf` typically signifies an uninitialized or invalid structure,
  while `Inf` often represents a null or discarded possibility during folding.
- `desc::String`: A human-readable description of the structural element. Common
  descriptions include:
  - `"HAIRPIN:<seq>"`: A hairpin loop closed by specific base pairs.
  - `"STACK:<bp1>/<bp2>"`: A stacked base pair.
  - `"BULGE:<size>"`: A bulge loop of a specific size.
  - `"INTERIOR_LOOP:<size1>/<size2>"`: An interior loop with specified sizes on each side.
  - `"BIFURCATION:<unpaired>n/<helices>h"`: A multibranch loop with a certain number of unpaired
    nucleotides and helices.
- `ij::Vector{Tuple{Int, Int}}`: A list of base-pair indices `(i, j)` involved in this
  structural element. For simple elements like a single hairpin or stack, this typically
  contains one tuple `(i, j)` indicating the 1-based indices of the paired nucleotides.
  For complex elements like multibranch loops, it can contain multiple tuples
  representing the various stems (helices) that constitute the junction.


# Examples

```jldoctest
julia> s1 = SeqFold.Structure() # Default, uninitialized
   0    0   -Inf

julia> s2 = SeqFold.Structure(-3.2, "HAIRPIN:CG/CG", [(2, 7)])
   2    7   -3.2  HAIRPIN:CG/CG

julia> s3 = SeqFold.Structure(5.1, "BIFURCATION:2n/3h", [(1, 10), (15, 20), (25, 30)])
   1   10    5.1  BIFURCATION:2n/3h

julia> Bool(s1) # Check if structure is valid (not -Inf)
false

julia> Bool(s2)
true

julia> s2 == SeqFold.Structure(-3.2, "HAIRPIN:CG/CG", [(2, 7)]) # Structures can be compared
true
```

# See also

[`fold`](@ref), [`dg`](@ref), [`dot_bracket`](@ref)
"""
mutable struct Structure
    e::Float64
    desc::String
    ij::Vector{Tuple{Int, Int}}

    function Structure(e::Float64=-Inf, desc::String="", ij::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[])
        return new(e, desc, ij)
    end
end


Base.:(==)(a::Structure, b::Structure) = a.e == b.e && a.ij == b.ij
Base.Bool(st::Structure) = !isinf(st.e) && st.e != -Inf
with_ij(st::Structure, ij::Vector{Tuple{Int, Int}}) = Structure(st.e, st.desc, ij)

function Base.show(io::IO, st::Structure)
    i, j = isempty(st.ij) ? (0, 0) : first(st.ij)
    @printf(io, "%4d %4d %6.1f  %-15s", i, j, st.e, st.desc)
end


const STRUCT_DEFAULT = Structure(-Inf)
const STRUCT_NULL = Structure(Inf)
const Structs = Vector{Vector{Structure}}


"""
    MeltingConditions(seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)

Represents buffer conditions for melting temperature calculation.

# Fields and Units
- `seq1_conc`, `seq2_conc`: Sequence concentrations in **nM** (nanomolar)
- `Na`, `K`, `Tris`, `Mg`, `dNTPs`: Buffer component concentrations in **mM** (millimolar)

## Examples
```jldoctest
julia> MeltingConditions(:pcr)
MeltingConditions (PCR preset)
  • NEB PCR buffer conditions for Taq DNA Polymerase
  • seq1 concentration: 250.0 nM (typical primer concentration)
  • seq2 concentration: 0.0 nM (asymmetric PCR)
  • Mg²⁺: 1.5 mM (optimal for Taq DNA Polymerase per NEB guidelines)
  • K⁺: 50.0 mM
  • Tris: 2.0 mM
  • dNTPs: 0.2 mM

julia> MeltingConditions(:std)
MeltingConditions (standard preset)
  • Standard hybridization buffer
  • seq1 concentration: 25.0 nM
  • seq2 concentration: 25.0 nM
  • Na⁺: 50.0 mM

julia> MeltingConditions(seq1_conc=10,seq2_conc=0,Na=0,K=0,Tris=10,Mg=0,dNTPs=0)
MeltingConditions (custom)
  • seq1 concentration: 10.0 nM
  • seq2 concentration: 0.0 nM
  • Tris: 10.0 mM

julia> MeltingConditions(10,0,0,0,0,5,0)
MeltingConditions (custom)
  • seq1 concentration: 10.0 nM
  • seq2 concentration: 0.0 nM
  • Mg²⁺: 5.0 mM (higher than NEB's recommended 1.5-2.0 mM range)

julia> MeltingConditions(seq1_conc=0,seq2_conc=0,Na=5,K=0,Tris=0,Mg=0,dNTPs=0)
ERROR: ArgumentError: DNA concentration is too low!
[...]

julia> MeltingConditions(seq1_conc=10,seq2_conc=10,Na=0,K=0,Tris=0,Mg=0,dNTPs=0)
ERROR: ArgumentError: No cations for salt correction (Mg=0.0, Na=0, K=0, Tris=0)
[...]
```
"""
struct MeltingConditions
    seq1_conc::Float64
    seq2_conc::Float64
    Na::Float64
    K::Float64
    Tris::Float64
    Mg::Float64
    dNTPs::Float64
    
    function MeltingConditions(
        seq1_conc::Real,
        seq2_conc::Real,
        Na::Real,
        K::Real,
        Tris::Real,
        Mg::Real,
        dNTPs::Real
    )
        if !all(>=(0.0), (seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs))
            throw(DomainError("Negative concentration detected!"))
        end
        if isapprox(seq1_conc+seq2_conc, 0.0, atol=1e-3)
            throw(ArgumentError("DNA concentration is too low!"))
        end
        if isapprox(Mg, 0.0, atol=1e-4) && isapprox((Na + K + Tris), 0.0, atol=1e-4)
            throw(ArgumentError("No cations for salt correction (Mg=$(Mg), Na=$(Na), K=$(K), Tris=$(Tris))"))
        end

        new(
            Float64(seq1_conc),
            Float64(seq2_conc),
            Float64(Na),
            Float64(K),
            Float64(Tris),
            Float64(Mg),
            Float64(dNTPs)
        )
    end
end

MeltingConditions(s::Symbol) = MeltingConditions(Val(s))
MeltingConditions(::Val{:pcr}) = PCR_CONDITIONS
MeltingConditions(::Val{:std}) = STD_CONDITIONS
MeltingConditions(t::Tuple{Vararg{Real,7}}) = MeltingConditions(t...)
MeltingConditions(; kwargs...) = MeltingConditions(NamedTuple(kwargs))
function MeltingConditions(nt::NamedTuple)
    MeltingConditions(nt.seq1_conc, nt.seq2_conc, nt.Na, nt.K, nt.Tris, nt.Mg, nt.dNTPs)
end

function MeltingConditions(base::MeltingConditions; kwargs...)
    fields = fieldnames(MeltingConditions)
    vals = [getfield(base, f) for f in fields]
    
    for (k, v) in pairs(kwargs)
        idx = findfirst(f -> string(f) == string(k), fields)
        if !isnothing(idx)
            vals[idx] = v
        else
            throw(ArgumentError("Unknown condition parameter: $k"))
        end
    end
    
    return MeltingConditions(vals...)
end

const PCR_CONDITIONS = MeltingConditions(
    seq1_conc = 250.0,
    seq2_conc = 0.0,
    Na = 0.0,
    K = 50.0,
    Tris = 2.0,
    Mg = 1.5,
    dNTPs = 0.2
)

const STD_CONDITIONS = MeltingConditions(
    seq1_conc = 25.0,
    seq2_conc = 25.0,
    Na = 50.0,
    K = 0.0,
    Tris = 0.0,
    Mg = 0.0,
    dNTPs = 0.0
)

function Base.show(io::IO, mc::MeltingConditions)
    if mc == PCR_CONDITIONS
        print(io, "MeltingConditions(:pcr)")
    elseif mc == STD_CONDITIONS
        print(io, "MeltingConditions(:std)")
    else
        params = []
        if mc.seq1_conc > 0
            push!(params, "seq1=$(round(mc.seq1_conc, digits=1))nM")
        end
        if mc.Mg > 0
            push!(params, "Mg=$(round(mc.Mg, digits=1))mM")
        end
        if mc.Na > 0
            push!(params, "Na=$(round(mc.Na, digits=1))mM")
        end
        
        if isempty(params)
            print(io, "MeltingConditions(custom)")
        else
            print(io, "MeltingConditions(", join(params, ","), ")")
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", mc::MeltingConditions)
    # Check if it matches a known preset
    if mc == PCR_CONDITIONS
        print(io, "MeltingConditions (PCR preset)\n")
        print(io, "  • NEB PCR buffer conditions for Taq DNA Polymerase\n")
        print(io, "  • seq1 concentration: 250.0 nM (typical primer concentration)\n")
        print(io, "  • seq2 concentration: 0.0 nM (asymmetric PCR)\n")
        print(io, "  • Mg²⁺: 1.5 mM (optimal for Taq DNA Polymerase per NEB guidelines)\n")
        print(io, "  • K⁺: 50.0 mM\n")
        print(io, "  • Tris: 2.0 mM\n")
        print(io, "  • dNTPs: 0.2 mM")
    elseif mc == STD_CONDITIONS
        print(io, "MeltingConditions (standard preset)\n")
        print(io, "  • Standard hybridization buffer\n")
        print(io, "  • seq1 concentration: 25.0 nM\n")
        print(io, "  • seq2 concentration: 25.0 nM\n")
        print(io, "  • Na⁺: 50.0 mM")
    else
        print(io, "MeltingConditions (custom)\n")
        print(io, "  • seq1 concentration: ", round(mc.seq1_conc, digits=1), " nM\n")
        print(io, "  • seq2 concentration: ", round(mc.seq2_conc, digits=1), " nM")
        
        # Only show non-zero cations with NEB context where appropriate
        if mc.Na > 0 
            print(io, "\n  • Na⁺: ", round(mc.Na, digits=1), " mM")
            if isapprox(mc.Na, 50.0, atol=0.1)
                print(io, " (standard hybridization buffer)")
            end
        end
        if mc.K > 0 
            print(io, "\n  • K⁺: ", round(mc.K, digits=1), " mM")
            if isapprox(mc.K, 50.0, atol=0.1)
                print(io, " (NEB Taq buffer component)")
            end
        end
        if mc.Tris > 0 
            print(io, "\n  • Tris: ", round(mc.Tris, digits=1), " mM")
            if isapprox(mc.Tris, 2.0, atol=0.1)
                print(io, " (NEB Taq buffer component)")
            end
        end
        if mc.Mg > 0 
            print(io, "\n  • Mg²⁺: ", round(mc.Mg, digits=1), " mM")
            if 1.5 <= mc.Mg <= 2.0
                print(io, " (optimal for Taq DNA Polymerase per NEB guidelines)")
            elseif mc.Mg > 2.0
                print(io, " (higher than NEB's recommended 1.5-2.0 mM range)")
            end
        end
        if mc.dNTPs > 0 
            print(io, "\n  • dNTPs: ", round(mc.dNTPs, digits=1), " mM")
            if isapprox(mc.dNTPs, 0.2, atol=0.05)
                print(io, " (NEB Taq buffer component)")
            end
        end
    end
end