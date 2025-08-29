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

mutable struct Struct
    e::Float64
    desc::String
    ij::Vector{Tuple{Int, Int}}

    function Struct(e::Float64=-Inf, desc::String="", ij::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[])
        return new(e, desc, ij)
    end
end


Base.:(==)(a::Struct, b::Struct) = a.e == b.e && a.ij == b.ij
Base.Bool(st::Struct) = !isinf(st.e) && st.e != -Inf
with_ij(st::Struct, ij::Vector{Tuple{Int, Int}}) = Struct(st.e, st.desc, ij)

function Base.show(io::IO, st::Struct)
    i, j = isempty(st.ij) ? (0, 0) : first(st.ij)
    @printf(io, "%4d %4d %6.1f  %-15s", i, j, st.e, st.desc)
end


const STRUCT_DEFAULT = Struct(-Inf)
const STRUCT_NULL = Struct(Inf)
const Structs = Vector{Vector{Struct}}