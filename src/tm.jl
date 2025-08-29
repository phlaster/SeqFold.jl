"""
    tm(seq1, seq2; pcr=true, concentrations=nothing) -> Float64
    tm(seq; pcr=true, concentrations=nothing) -> Float64

Calculate the melting temperature (Tm) of a DNA sequence using nearest-neighbor thermodynamics.

Compute the melting temperature (°C) of a DNA duplex. The first method calculates Tm for two sequences
of the same length, while the second method uses the complement of the input sequence as the second strand.

# Arguments
- `seq1::AbstractString`: DNA sequence
- `seq2::AbstractString`: Another DNA sequence, must be the same length as `seq1`
- `seq::AbstractString`: Single DNA sequence to be matched with its exact complement
- `pcr::Bool=true`: Whether to use PCR conditions (affects salt concentration assumptions)
- `concentrations`: should be either `nothing` (`pcr` flag is used in such case) or
    `(seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)::NTuple{7, Float64}`
    providing exact concentrations (`pcr` flag is ignored in such case).

# Examples
```jldoctest
julia> tm("CACACTTAGGGTAGCATCGA")
61.3

julia> tm("CACACTTAGGGTAGCATCGA", "ATCGCTTGAGGTAGCGTTGA")
-275.6

julia> tm("CACACTTAGGGTAGCATCGA", pcr=false)
50.7
```
# Implementation
The calculation uses nearest-neighbor thermodynamic parameters from related literature (follow links below),
accounting for initialization terms, nearest-neighbor pairs, and terminal mismatches when present.

# See also
[`DNA_NN`](@ref), [`DNA_INTERNAL_MM`](@ref), [`DNA_TERMINAL_MM`](@ref)

"""
function tm(seq1::AbstractString, seq2::AbstractString;
    pcr::Bool=true,
    concentrations::Union{Nothing, NTuple{7, Float64}}=nothing
    )::Float64
    seq1, seq2, seq_len, gc_count = _parse_input(seq1, seq2)
    gc = gc_count / seq_len

    dh, ds = DNA_ENERGIES.NN["init"]
    
    init = string(first(seq1), last(seq1))
    init_at = count(c -> c in "AT", init)
    init_gc = count(c -> c in "GC", init)
    
    init_at_h, init_at_s = DNA_ENERGIES.NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_ENERGIES.NN["init_G/C"]
    dh += init_at * init_at_h + init_gc * init_gc_h
    ds += init_at * init_at_s + init_gc * init_gc_s
    
    @inbounds for i in 1:seq_len-1
        pair = string(seq1[i], seq1[i+1], '/', seq2[i], seq2[i+1])

        pair_dh, pair_ds = if (i == 1 || i == seq_len-1) && haskey(DNA_ENERGIES.TERMINAL_MM, pair)
            DNA_ENERGIES.TERMINAL_MM[pair]
        elseif haskey(DNA_ENERGIES.NN, pair)
            DNA_ENERGIES.NN[pair]
        elseif haskey(DNA_ENERGIES.INTERNAL_MM, pair)
            DNA_ENERGIES.INTERNAL_MM[pair]
        else
            0.0, 0.0
        end
        
        dh += pair_dh
        ds += pair_ds
    end
    
    if isnothing(concentrations)
        return _calc_tm(dh, ds, pcr, gc, seq_len)
    else
        return _calc_tm(dh, ds, gc, seq_len, concentrations)
    end
end

tm(seq::AbstractString; pcr::Bool=true, concentrations=nothing) = tm(seq, complement(seq); pcr=pcr, concentrations=concentrations)

function _parse_input(s1::AbstractString, s2::AbstractString)
    seq1, n1, gc_count1 = _parse_input(s1)
    seq2, n2, _ = _parse_input(s2)

    @assert n1 == n2 "Length mismatch between sequences: $n1 vs $n2"

    return seq1, seq2, n1, gc_count1
end

function _parse_input(seq::AbstractString)
    n = length(seq)
    @assert n >= 2 "Sequences lengths have to be at least 2bp"
    
    seq = uppercase(seq)
    gc_count = 0

    for c in seq
        @assert in(c, "ACGT") "Incorrect character '$c' in DNA sequence"
        gc_count += Int(c == 'G' || c == 'C')
    end

    return seq, n, gc_count
end

"""
    tm_cache(seq; pcr=true) -> Matrix{Float64}
    tm_cache(seq1, seq2; pcr=true) -> Matrix{Float64}

Compute a matrix of melting temperatures for all possible subsequences of a DNA sequence pair.

The resulting matrix has element (`i`, `j`) representing the melting temperature of the DNA subsequence
from position i to position j (1-based indexing). This cache enables efficient Tm calculations for
various subsequences without recalculating thermodynamic parameters.

# Arguments
- `seq1::AbstractString`: DNA sequence
- `seq2::AbstractString`: Another DNA sequence, must be the same length as `seq1`
- `seq::AbstractString`: Single DNA sequence to be matched with its exact complement
- `pcr::Bool=true`: Whether to use PCR conditions (affects salt concentration assumptions)

# Returns
A `Matrix{Float64}` where element (`i`, `j`) contains the melting temperature (in °C) of the subsequence
from position `i` to position `j`, inclusive. Elements where `j < i` contain `Inf` as they represent
invalid ranges, and single-nucleotide subsequences also have `Inf` as they don't have meaningful Tm values.

# Errors
Throws an `ErrorException` if a DNA pair (e.g., "AT/CC") is not found in the nearest-neighbor energy tables.
The energy tables contain only standard Watson-Crick base pairs as documented related papers.
Non-standard pairings or invalid sequences will cause this error. See links below for source papers.

# Examples
```jldoctest
julia> tm_cache("ATGCC")
5×5 Matrix{Float64}:
 Inf  -66.2  -20.9    2.1   18.4
 Inf   Inf   -45.2  -11.9    9.7
 Inf   Inf    Inf   -43.1   -9.4
 Inf   Inf    Inf    Inf   -52.7
 Inf   Inf    Inf    Inf    Inf

julia> tm_cache("AAGC", "TTCG")
4×4 Matrix{Float64}:
 Inf  -67.8  -21.7    5.6
 Inf   Inf   -48.8   -8.7
 Inf   Inf    Inf   -34.6
 Inf   Inf    Inf    Inf

julia> tm_cache("AAGC", "TTCG"; pcr=false)
4×4 Matrix{Float64}:
 Inf   0.000202691   0.000202305   0.000202048
 Inf  Inf            0.000202048   0.000201791
 Inf  Inf           Inf            0.000201406
 Inf  Inf           Inf           Inf

julia> tm_cache("AAGC", "TAAG")
ERROR: Pair AG/AA not found in energy tables
[...]
```


# Implementation
The function uses dynamic programming to build a cache of Tm values for all subsequences.
The algorithm has O(n²) time complexity and O(n²) space complexity, where n is the sequence length.
This approach avoids redundant calculations when multiple Tm values for different subsequences are needed.

# See also
[`tm`](@ref), [`gc_cache`](@ref), [`DNA_NN`](@ref), [`DNA_INTERNAL_MM`](@ref)
"""
function tm_cache(seq1::AbstractString, seq2::AbstractString; pcr::Bool=true)::Matrix{Float64}
    seq1, seq2, n, _ = _parse_input(seq1, seq2)
    
    arr_gc = gc_cache(seq1)
    arr_dh = zeros(n, n)
    arr_ds = zeros(n, n)
    arr_tm = fill(Inf, n, n)
    
    @inbounds for i in 1:n
        if i == n
            arr_dh[i, i] = arr_dh[i-1, i-1]
            arr_ds[i, i] = arr_ds[i-1, i-1]
            continue
        end
        
        pair = string(seq1[i], seq1[i+1], '/', seq2[i], seq2[i+1])
        dh, ds = if haskey(DNA_ENERGIES.NN, pair)
            DNA_ENERGIES.NN[pair]
        elseif haskey(DNA_ENERGIES.INTERNAL_MM, pair)
            DNA_ENERGIES.INTERNAL_MM[pair]
        else
            error("Pair $pair not found in energy tables")
        end

        arr_dh[i, i] = dh
        arr_ds[i, i] = ds
    end
    
    for i in 1:n
        for j in i+1:n
            arr_dh[i, j] = arr_dh[i, j-1] + arr_dh[j, j]
            arr_ds[i, j] = arr_ds[i, j-1] + arr_ds[j, j]
            len = j - i + 1
            arr_tm[i, j] = _calc_tm(arr_dh[i, j], arr_ds[i, j], pcr, arr_gc[i, j], len)
        end
    end
    
    return arr_tm
end

tm_cache(seq::AbstractString; pcr::Bool=true)::Matrix{Float64} = tm_cache(seq, complement(seq); pcr=pcr)

"""
    gc_cache(seq) -> Matrix{Float64}

Compute a matrix of GC scores for all possible subsequences of a DNA sequence.

The resulting matrix has element `[i, j]` representing the GC score
of the DNA subsequence from position `i` to position `j`. This cache enables efficient
GC ratio calculations for various subsequences without redundant computations.

# Arguments
- `seq::AbstractString`: The DNA sequence to analyze

# Returns
A `Matrix{Float64}` where element `[i, j]` contains the GC ratio of the subsequence from position `i` to `j`,
rounded to one decimal place. Elements where `j < i` contain `Inf` as they represent invalid ranges.

# Examples
```jldoctest
julia> gc_cache("GGAA")
4×4 Matrix{Float64}:
  1.0   1.0   0.7  0.5
 Inf    1.0   0.5  0.3
 Inf   Inf    0.0  0.0
 Inf   Inf   Inf   0.0

julia> gc_cache("GCTA")
4×4 Matrix{Float64}:
  1.0   1.0   0.7  0.5
 Inf    1.0   0.5  0.3
 Inf   Inf    0.0  0.0
 Inf   Inf   Inf   0.0

julia> gc_cache("ATA")
3×3 Matrix{Float64}:
  0.0   0.0  0.0
 Inf    0.0  0.0
 Inf   Inf   0.0

julia> gc_cache("GGTT") == gc_cache("CCAA")
true
```

# Implementation
The algorithm has O(n²) time complexity and O(n²) space complexity, where `n` is the sequence length.
This approach avoids redundant calculations when multiple GC ratios for different subsequences are needed.

# See also
[`tm_cache`](@ref), [`tm`](@ref)
"""
function gc_cache(seq::AbstractString)::Matrix{Float64}
    seq, n, _ = _parse_input(seq)
    arr_gc = fill(Inf, n, n)
    
    @inbounds for i in 1:n
        if i == n
            arr_gc[i, i] = arr_gc[i-1, i-1]
            continue
        end
        arr_gc[i, i] = Float64(seq[i] in "GC")
        if i == n-1 && iszero(arr_gc[i, i])
            arr_gc[i, i] = Float64(seq[i+1] in "GC")
        end
    end
    
    @inbounds for i in 1:n
        @simd for j in i+1:n
            arr_gc[i, j] = arr_gc[i, j-1] + arr_gc[j, j]
        end
    end
    
    @inbounds for i in 1:n
        @simd for j in i:n
            len = j - i + 1
            arr_gc[i, j] = round(arr_gc[i, j] / len, digits=1)
        end
    end
    
    return arr_gc
end

function _calc_tm(dh::Float64, ds::Float64, pcr::Bool, gc::Float64, seq_len::Int)::Float64
    seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs = pcr ?
    (250.0,  0.0,  0.0, 50.0, 2.0, 1.5, 0.2) :
    ( 25.0, 25.0, 50.0,  0.0, 0.0, 0.0, 0.0)
    return _calc_tm(dh, ds, gc, seq_len, seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)
end

function _calc_tm(
    dh::Float64, ds::Float64, gc::Float64, seq_len::Int,
    (seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)::NTuple{7, Float64}
    )::Float64
    return _calc_tm(dh, ds, gc, seq_len, seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)
end

function _calc_tm(dh::Float64, ds::Float64, gc::Float64, seq_len::Int,
    seq1_conc::Real,
    seq2_conc::Real,
    Na::Real,
    K::Real,
    Tris::Real,
    Mg::Real,
    dNTPs::Real,
    )::Float64

    # salt correction for deltaS
    # copied-pasted from Bio.SeqUtils' use of a decision tree by:
    # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
    Mon = Na + K + Tris / 2.0
    mg = Mg * 1e-3 # Lowercase ions (mg, mon, dntps) are molar
    mon = Mon * 1e-3
    
    # coefficients to a multi-variate from the paper
    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31
    
    if dNTPs > 0.0
        dntps = dNTPs * 1e-3
        ka = 3e4 # Dissociation constant for Mg:dNTP
        # Free Mg2+ calculation:
        mg = (-(ka * dntps - ka * mg + 1.0) + sqrt((ka * dntps - ka * mg + 1.0)^2 + 4.0 * ka * mg)) / (2.0 * ka)
    end
    
    if isapprox(mg, 0.0, atol=1e-7)
        if isapprox(mon, 0.0, atol=1e-7)
            error("No cations for salt correction")
        end
        corr = (4.29 * gc - 3.95) * 1e-5 * log(mon) + 9.4e-6 * log(mon)^2
    else
        if isapprox(mon, 0.0, atol=1e-7)
            # Use defaults as if R is large
            corr = (
                a
                + b * log(mg)
                + gc * (c + d * log(mg))
                + (inv(2.0 * (seq_len - 1))) * (e + f * log(mg) + g * log(mg) ^ 2)
            ) * 1e-5
        else
            R = sqrt(mg) / mon
            if R < 0.22
                corr = (4.29 * gc - 3.95) * 1e-5 * log(mon) + 9.4e-6 * log(mon) ^ 2
            else
                if R < 6.0
                    a = 3.92 * (0.843 - 0.352 * sqrt(mon) * log(mon))
                    d = 1.42 * (1.279 - 4.03e-3 * log(mon) - 8.03e-3 * log(mon) ^ 2)
                    g = 8.31 * (0.486 - 0.258 * log(mon) + 5.25e-3 * log(mon) ^ 3)
                end
                corr = (
                    a
                    + b * log(mg)
                    + gc * (c + d * log(mg))
                    + (inv(2.0 * (seq_len - 1))) * (e + f * log(mg) + g * log(mg) ^ 2)
                ) * 1e-5
            end
        end
    end
    
    # tm with concentration consideration
    k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9
    R_const = 1.9872
    est = (dh * 1000.0) / (ds + R_const * log(k)) - 273.15

    # add in salt correction
    est = inv(inv(est + 273.15) + corr) - 273.1
    
    return round(max(est, -273.2), digits=1)
end

_comp_dna(b::UInt8)::UInt8 =
    b == UInt8('A') ? UInt8('T') :
    b == UInt8('a') ? UInt8('T') :
    b == UInt8('T') ? UInt8('A') :
    b == UInt8('t') ? UInt8('A') :
    b == UInt8('C') ? UInt8('G') :
    b == UInt8('c') ? UInt8('G') :
    b == UInt8('G') ? UInt8('C') :
    b == UInt8('g') ? UInt8('C') :
    UInt8('.')

const DNA_COMP_TABLE = collect(ntuple(i->_comp_dna(UInt8(i-1)), Val(256)))

@inline function _comp_bytes!(out::Vector{UInt8}, seq::AbstractVector{UInt8})
    n = length(seq)
    @inbounds @simd for i in 1:n
        out[i] = DNA_COMP_TABLE[ seq[i] + 1 ]
    end
    return out
end

function complement(seq::AbstractVector{UInt8})
    out = Vector{UInt8}(undef, length(seq))
    return _comp_bytes!(out, seq)
end

function complement(seq::AbstractString)
    out_bytes = complement(codeunits(seq))
    return String(out_bytes)
end