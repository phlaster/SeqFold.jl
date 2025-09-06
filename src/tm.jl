"""
    tm(seq1, seq2; conditions=:pcr, kwargs...) -> Float64
    tm(seq; conditions=:pcr, kwargs...)  -> Float64

Calculate melting temperature (Tm °C) for DNA duplex formation using nearest-neighbor thermodynamics.

# Arguments
- `seq1::AbstractString`: DNA sequence;
- `seq2::AbstractString`: Another DNA sequence, must be the same length as `seq1`;
- `seq::AbstractString`: Single DNA sequence to be matched with its exact complement;
- `conditions`: Buffer conditions specification (see below);
- `kwargs...`: Additional parameters to override preset conditions in place.

# Conditions Specification
The `conditions` parameter can be:
- `:pcr` (default) or `:std`: Use preset conditions;
- `MeltingConditions` object: Custom conditions (for more info see: [`MeltingConditions`](@ref));
- `NTuple{7, Float64}`: `(seq1_conc, seq2_conc, Na, K, Tris, Mg, dNTPs)`;
- `NamedTuple` with condition fields.

# Examples
```jldoctest
julia> tm("GGGGGG")
15.5

julia> tm("GGGGGG", conditions=:pcr)
15.5

julia> tm("GGGGGG", Mg=10)
23.0

julia> tm("GGGGGG", conditions=:std)
2.6

julia> tm("GGGGGG", conditions=:std, Na=100)
5.1

julia> tm("GGGGGG", conditions=MeltingConditions(150, 150, 20, 0, 0, 10, 0))
16.8

julia> tm("ACCCCC", "GGGGGG")
6.2
```
# Implementation
The calculation uses nearest-neighbor thermodynamic parameters from related literature (follow links below to see sources),
accounting for initialization terms, nearest-neighbor pairs, and terminal mismatches when present.

# See also
[`MeltingConditions`](@ref)
"""
function tm(seq1::AbstractString, seq2::AbstractString; conditions=:pcr, kwargs...)
    base_cond = MeltingConditions(conditions)
    cond = isempty(kwargs) ? base_cond : MeltingConditions(base_cond; kwargs...)
    
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
            (0.0, 0.0)
        end
        
        dh += pair_dh
        ds += pair_ds
    end
    
    return _calc_tm(dh, ds, gc, seq_len, cond)
end

function tm(seq::AbstractString; conditions=:pcr, kwargs...)
    tm(seq, complement(seq); conditions=conditions, kwargs...)
end

"""
    tm_cache(seq1, seq2; conditions=:pcr, kwargs...) -> Matrix{Float64}
    tm_cache(seq; conditions=:pcr, kwargs...) -> Matrix{Float64}

Compute a matrix of melting temperatures for all possible subsequences of a DNA sequence pair.

# Arguments
- `seq1::AbstractString`: DNA sequence;
- `seq2::AbstractString`: Another DNA sequence, must be the same length as `seq1`;
- `seq::AbstractString`: Single DNA sequence to be matched with its exact complement;
- `conditions`: Buffer conditions specification (see `tm` docstrings);
- `kwargs...`: Additional parameters to override preset conditions in place.

# Returns
A `Matrix{Float64}` where element (`i`, `j`) contains the melting temperature (in °C) of the subsequence
from position `i` to position `j`, inclusive. Elements where `j < i` contain `Inf` as they represent
invalid ranges, and single-nucleotide subsequences also have `Inf` as they don't have meaningful Tm values.

# Examples
```jldoctest
julia> SeqFold.tm_cache("ATCC")
4×4 Matrix{Float64}:
 Inf  -212.6   -95.3   -48.6
 Inf    Inf   -161.6   -82.7
 Inf    Inf     Inf   -135.5
 Inf    Inf     Inf     Inf

julia> SeqFold.tm_cache("AAGC", "TTCG")
4×4 Matrix{Float64}:
 Inf  -204.8   -94.6   -40.3
 Inf    Inf   -166.7   -72.9
 Inf    Inf     Inf   -116.7
 Inf    Inf     Inf     Inf

julia> SeqFold.tm_cache("AAGC", "TTCG"; conditions=:std)
4×4 Matrix{Float64}:
 Inf  -213.1  -109.3   -55.6
 Inf    Inf   -177.9   -88.1
 Inf    Inf     Inf   -129.8
 Inf    Inf     Inf     Inf
```

# Implementation
The function uses dynamic programming to build a cache of Tm values for all subsequences.
The algorithm has `O(n²)` time and space complexity, where `n` is the sequence length.
This approach avoids redundant calculations when multiple Tm values for different subsequences are needed.

# See also
[`tm`](@ref), [`SeqFold.gc_cache`](@ref)
"""
function tm_cache(seq1::AbstractString, seq2::AbstractString; conditions=:pcr, kwargs...)::Matrix{Float64}
    base_cond = MeltingConditions(conditions)
    cond = isempty(kwargs) ? base_cond : MeltingConditions(base_cond; kwargs...)
    seq1, seq2, n, _ = _parse_input(seq1, seq2)
    
    base_dh, base_ds = DNA_ENERGIES.NN["init"]
    init_at_h, init_at_s = DNA_ENERGIES.NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_ENERGIES.NN["init_G/C"]
    
    term_pair_h = zeros(Float64, n-1)
    term_pair_s = zeros(Float64, n-1)
    non_term_pair_h = zeros(Float64, n-1)
    non_term_pair_s = zeros(Float64, n-1)
    
    @inbounds for k in 1:(n-1)
        pair = string(seq1[k], seq1[k+1], '/', seq2[k], seq2[k+1])
        
        if haskey(DNA_ENERGIES.NN, pair)
            non_term_pair_h[k], non_term_pair_s[k] = DNA_ENERGIES.NN[pair]
        elseif haskey(DNA_ENERGIES.INTERNAL_MM, pair)
            non_term_pair_h[k], non_term_pair_s[k] = DNA_ENERGIES.INTERNAL_MM[pair]
        else
            non_term_pair_h[k], non_term_pair_s[k] = 0.0, 0.0
        end
        
        if haskey(DNA_ENERGIES.TERMINAL_MM, pair)
            term_pair_h[k], term_pair_s[k] = DNA_ENERGIES.TERMINAL_MM[pair]
        else
            term_pair_h[k], term_pair_s[k] = non_term_pair_h[k], non_term_pair_s[k]
        end
    end
    
    prefix_non_term_h = zeros(Float64, n)
    prefix_non_term_s = zeros(Float64, n)
    
    @inbounds for k in 2:n
        prefix_non_term_h[k] = prefix_non_term_h[k-1] + non_term_pair_h[k-1]
        prefix_non_term_s[k] = prefix_non_term_s[k-1] + non_term_pair_s[k-1]
    end
    
    gc_matrix = gc_cache(seq1)
    
    cache = fill(Inf, n, n)
    
    @inbounds for L in 2:n
        for i in 1:(n - L + 1)
            j = i + L - 1
            
            init = string(seq1[i], seq1[j])
            init_at = count(c -> c in "AT", init)
            init_gc = count(c -> c in "GC", init)
            dh = base_dh + init_at * init_at_h + init_gc * init_gc_h
            ds = base_ds + init_at * init_at_s + init_gc * init_gc_s
            
            dh += prefix_non_term_h[j] - prefix_non_term_h[i]
            ds += prefix_non_term_s[j] - prefix_non_term_s[i]
            
            dh += term_pair_h[i] - non_term_pair_h[i]
            ds += term_pair_s[i] - non_term_pair_s[i]
            
            if L > 2
                dh += term_pair_h[j-1] - non_term_pair_h[j-1]
                ds += term_pair_s[j-1] - non_term_pair_s[j-1]
            end
            
            gc = gc_matrix[i,j]
            cache[i, j] = _calc_tm(dh, ds, gc, L, cond)
        end
    end
    
    return cache
end

function tm_cache(seq::AbstractString; conditions=:pcr, kwargs...)::Matrix{Float64}
    tm_cache(seq, complement(seq); conditions=conditions, kwargs...)
end

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
julia> SeqFold.gc_cache("GGAA")
4×4 Matrix{Float64}:
  1.0   1.0   0.666667  0.5
 Inf    1.0   0.5       0.333333
 Inf   Inf    0.0       0.0
 Inf   Inf   Inf        0.0

julia> SeqFold.gc_cache("GAAA")
4×4 Matrix{Float64}:
  1.0   0.5   0.333333  0.25
 Inf    0.0   0.0       0.0
 Inf   Inf    0.0       0.0
 Inf   Inf   Inf        0.0

julia> SeqFold.gc_cache("ATA")
3×3 Matrix{Float64}:
  0.0   0.0  0.0
 Inf    0.0  0.0
 Inf   Inf   0.0

julia> SeqFold.gc_cache("GGTT") == SeqFold.gc_cache("CCAA")
true
```

# See also
[`SeqFold.tm_cache`](@ref), [`tm`](@ref)
"""
function gc_cache(seq::AbstractString)::Matrix{Float64}
    n = length(seq)
    prefix = zeros(Int, n+1)
    @inbounds @simd for i in 1:n
        prefix[i+1] = prefix[i] + (seq[i] in "GC" ? 1 : 0)
    end
    
    M = fill(Inf, n, n)
    
    @inbounds @simd for i in 1:n
        for j in i:n
            len = j - i + 1
            gc_count = prefix[j+1] - prefix[i]
            M[i, j] = gc_count / len
        end
    end
    
    return M
end

function _calc_tm(
    dh::Float64, 
    ds::Float64, 
    gc::Float64, 
    seq_len::Int,
    cond::MeltingConditions
    )::Float64

    # salt correction for deltaS
    # copied-pasted from Bio.SeqUtils' use of a decision tree by:
    # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
    Mon = cond.Na + cond.K + cond.Tris / 2.0
    mg = cond.Mg * 1e-3 # Lowercase ions (mg, mon, dntps) are molar
    mon = Mon * 1e-3
    
    # coefficients to a multi-variate from the paper
    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31
    
    if cond.dNTPs > 0.0
        dntps = cond.dNTPs * 1e-3
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
    k = abs((cond.seq1_conc - (cond.seq2_conc / 2.0)) * 1e-9) + 1e-15
    R_const = 1.9872
    inv_est = (ds + R_const * log(k)) / (dh * 1000.0)

    # add in salt correction
    est = inv(inv_est + corr) - 273.15
    
    return round(max(est, -273.2), digits=1)
end
