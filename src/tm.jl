const IUPAC_MAP = Dict(
    'A'=>"A",  'C'=>"C",  'G'=>"G",  'T'=>"T",
    'R'=>"AG", 'Y'=>"CT", 'S'=>"CG", 'W'=>"AT",
    'K'=>"GT", 'M'=>"AC", 'B'=>"CGT",'D'=>"AGT",
    'H'=>"ACT",'V'=>"ACG",'N'=>"ACGT"
)
const IUPAC_COUNTS = Dict(
    'A' => 1, 'C' => 1, 'G' => 1, 'T' => 1,
    'M' => 2, 'R' => 2, 'W' => 2, 'S' => 2, 'Y' => 2, 'K' => 2,
    'V' => 3, 'H' => 3, 'D' => 3, 'B' => 3, 'N' => 4
)

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
    tm_deg(seq; conditions=:pcr, high_deg_warn=10^5, conf_level=0.9, kwargs...) -> NamedTuple

Calculate statistics of melting temperature (Tm °C) for a degenerate DNA sequence by enumerating all possible 
non-degenerate variants and computing their melting temperatures using nearest-neighbor thermodynamics.

This function handles IUPAC ambiguity codes by generating all possible non-degenerate 
sequences, computing [`tm`](@ref) for each, and returning both the average temperature and a confidence interval 
that contains `conf_level` proportion of variants. A warning is issued when the number of 
variants exceeds `high_deg_warn` as computation may become slow for highly degenerate sequences.

# Arguments
- `seq::AbstractString`: Degenerate DNA sequence (with IUPAC ambiguity codes);
- `conditions`: Buffer conditions specification (see [`tm`](@ref) for detailed review);
- `high_deg_warn`: Threshold for warning about high degeneracy (number of variants). Set non-positive or `false` to disable;
- `conf_level::Float64`: Confidence level `(0,1]` for the interval. For example, 0.95 means 95% of variants fall within the interval;
- `kwargs...`: Additional parameters to override preset conditions in place.

# Returns
A `NamedTuple` with fields:
- `mean::Float64`: Average melting temperature across all variants (rounded to 1 decimal place)
- `conf::Tuple{Float64, Float64}`: Confidence interval bounds `(lower, upper)` containing `conf_level` proportion of variants

# Examples
```jldoctest
julia> tm_deg("GGGGGG") # non-degenerate sequence
(mean = 15.5, conf = (15.5, 15.5))

julia> result = tm_deg("GGGGGW") # W for A or T
(mean = 10.5, conf = (10.0, 11.0))

julia> T_a, T_t = tm("GGGGGA"), tm("GGGGGT")
(10.0, 11.0)

julia> isapprox(result.mean, (T_a+T_t)/2, atol=0.1)
true

julia> result.conf == (min(T_a, T_t), max(T_a, T_t)) # for binary degeneracy, interval spans all values
true

julia> tm_deg("GGGGGW", conditions=:std) # standard conditions
(mean = -3.4, conf = (-3.9, -2.9))

julia> tm_deg("GGGGGW", Mg=5) # adjust conditions in-place
(mean = 16.2, conf = (15.7, 16.8))

julia> tm_deg("NNNNNNN", conf_level=0.5) # 50% confidence interval (middle 50% of variants)
(mean = 11.3, conf = (4.8, 17.7))

julia> tm_deg("NNNNNNN", conf_level=0.9) # 90% confidence interval
(mean = 11.3, conf = (-4.4, 27.1))

julia> tm_deg("NNNNNNNNN");
┌ Warning: High degeneracy: 9 positions → 262144 variants. Computation may be slow.
└ @ SeqFold ~/work/SeqFold.jl/SeqFold.jl/src/tm.jl:210

julia> tm_deg("NNNNNNNNN", high_deg_warn=1e6) # increase warning threshold
(mean = 28.4, conf = (14.8, 42.2))

julia> tm_deg("NNNNNNNNN", high_deg_warn=false) # disable warning entirely
(mean = 28.4, conf = (14.8, 42.2))
```

# Note!
Recommended for sequences with ≤5 degenerate positions; for higher degeneracy, computation becomes slow.

# See also
[`tm`](@ref), [`MeltingConditions`](@ref)
"""
function tm_deg(seq::AbstractString; 
                conditions=:pcr, 
                high_deg_warn=10^5,
                conf_level=0.9,
                kwargs...)::@NamedTuple{mean::Float64, conf::NTuple{2, Float64}}
    
    if !(0 < conf_level <= 1)
        throw(ArgumentError("conf_level must be between 0 and 1"))
    end
    
    base_cond = MeltingConditions(conditions)
    cond = isempty(kwargs) ? base_cond : MeltingConditions(base_cond; kwargs...)

    seq = uppercase(seq)
    allowed_chars = Set("ACGTMRWSYKVHDBN")
    seq_chars = Set(seq)
    if !issubset(seq_chars, allowed_chars)
        throw(ArgumentError("Input sequence contains unallowed characters: $(join(collect(setdiff(seq_chars, allowed_chars)), ", "))"))
    end

    n_degenerate = 0
    n_nondeg = BigInt(1)
    
    @simd for char in seq
        if char in "MRWSYKVHDBN"
            n_degenerate += 1
        end
        n_nondeg *= IUPAC_COUNTS[char]
    end

    if n_nondeg == 1
        t = tm(seq; conditions=cond)
        return (mean=t, conf=(t, t))
    end
    
    if n_nondeg > high_deg_warn && high_deg_warn > 0
        @warn "High degeneracy: $n_degenerate positions → $n_nondeg variants. Computation may be slow."
    end

    seq_len = length(seq)
    options = Vector{String}(undef, seq_len)
    lens = Vector{Int}(undef, seq_len)
    for (i, c) in enumerate(seq)
        options[i] = IUPAC_MAP[c]
        lens[i] = length(options[i])
    end

    indices = ones(Int, seq_len)
    buffer = Vector{Char}(undef, seq_len)
    
    tm_values = Vector{Float64}(undef, n_nondeg)

    @inbounds for i in 1:n_nondeg
        for j in 1:seq_len
            buffer[j] = options[j][indices[j]]
        end
        tm_values[i] = tm(String(buffer); conditions=cond)
        
        pos = seq_len
        while pos > 0
            indices[pos] += 1
            if indices[pos] <= lens[pos]
                break
            end
            indices[pos] = 1
            pos -= 1
        end
    end

    mean_tm = Float64(sum(tm_values) / n_nondeg)
    
    low_percentile = (1 - conf_level) / 2
    high_percentile = 1 - low_percentile
    
    sort!(tm_values)
    
    low_idx = max(1, min(n_nondeg, round(Int, low_percentile * (n_nondeg-1)) + 1))
    high_idx = max(1, min(n_nondeg, round(Int, high_percentile * (n_nondeg-1)) + 1))
    
    conf_interval = (
        round(tm_values[low_idx], digits=1), 
        round(tm_values[high_idx], digits=1)
    )
    
    return (
        mean = round(mean_tm, digits=1),
        conf = conf_interval
    )
end

"""
    SeqFold.tm_cache(seq1, seq2; conditions=:pcr, kwargs...) -> Matrix{Float64}
    SeqFold.tm_cache(seq; conditions=:pcr, kwargs...) -> Matrix{Float64}

Compute a matrix of melting temperatures for all possible subsequences of a DNA sequence pair.

# Arguments
- `seq1::AbstractString`: DNA sequence;
- `seq2::AbstractString`: Another DNA sequence, must be the same length as `seq1`;
- `seq::AbstractString`: Single DNA sequence to be matched with its exact complement;
- `conditions`: Buffer conditions specification (see `tm` docstrings);
- `kwargs...`: Additional parameters to override preset conditions in place.

# Returns
A `Matrix{Float64}` where element `[i, j]` contains the melting temperature (in °C) of the subsequence
from position `i` to position `j`, inclusive. Elements where `j < i` contain `NaN` as they represent
invalid ranges, and single-nucleotide subsequences also have `NaN` as they don't have meaningful Tm values.

# Examples
```jldoctest
julia> SeqFold.tm_cache("ATCC")
4×4 Matrix{Float64}:
 NaN  -212.6   -95.3   -48.6
 NaN   NaN    -161.6   -82.7
 NaN   NaN     NaN    -135.5
 NaN   NaN     NaN     NaN

julia> SeqFold.tm_cache("AAGC", "TTCG")
4×4 Matrix{Float64}:
 NaN  -204.8   -94.6   -40.3
 NaN   NaN    -166.7   -72.9
 NaN   NaN     NaN    -116.7
 NaN   NaN     NaN     NaN

julia> SeqFold.tm_cache("AAGC", "TTCG"; conditions=:std)
4×4 Matrix{Float64}:
 NaN  -213.1  -109.3   -55.6
 NaN   NaN    -177.9   -88.1
 NaN   NaN     NaN    -129.8
 NaN   NaN     NaN     NaN
```

# Implementation
The function uses dynamic programming to build a cache of Tm values for all subsequences.
The algorithm has `O(n²)` time and space complexity, where `n` is the sequence length.
This approach avoids redundant calculations when multiple Tm values for different subsequences are needed.

# See also
[`tm`](@ref), [`SeqFold.gc_cache`](@ref), [`SeqFold.dg_cache`](@ref)
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
    
    cache = fill(NaN64, n, n)
    
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
    SeqFold.gc_cache(seq) -> Matrix{Float64}

Compute a matrix of GC scores for all possible subsequences of a DNA sequence.

The resulting matrix has element `[i, j]` representing the GC score
of the DNA subsequence from position `i` to position `j`. This cache enables efficient
GC ratio calculations for various subsequences without redundant computations.

# Arguments
- `seq::AbstractString`: The DNA sequence to analyze

# Returns
A `Matrix{Float64}` where element `[i, j]` contains the GC ratio of the subsequence from position `i` to `j`.
Elements where `j < i` contain `NaN` as they represent invalid ranges.

# Examples
```jldoctest
julia> SeqFold.gc_cache("GGAA")
4×4 Matrix{Float64}:
   1.0    1.0    0.666667  0.5
 NaN      1.0    0.5       0.333333
 NaN    NaN      0.0       0.0
 NaN    NaN    NaN         0.0

julia> SeqFold.gc_cache("GAAA")
4×4 Matrix{Float64}:
   1.0    0.5    0.333333  0.25
 NaN      0.0    0.0       0.0
 NaN    NaN      0.0       0.0
 NaN    NaN    NaN         0.0

julia> SeqFold.gc_cache("ATA")
3×3 Matrix{Float64}:
   0.0    0.0  0.0
 NaN      0.0  0.0
 NaN    NaN    0.0
```

# See also
[`SeqFold.tm_cache`](@ref), [`tm`](@ref), [`SeqFold.dg_cache`](@ref)
"""
function gc_cache(seq::AbstractString)::Matrix{Float64}
    n = length(seq)
    prefix = zeros(Int, n+1)
    @inbounds @simd for i in 1:n
        prefix[i+1] = prefix[i] + (seq[i] in "GC" ? 1 : 0)
    end
    
    M = fill(NaN64, n, n)
    
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
