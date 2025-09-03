"""
    fold(seq; temp = 37.0) -> Vector{Structure}

Predict the minimum free energy secondary structure of a nucleic acid sequence using a
dynamic programming algorithm based on the Zuker and Stiegler (1981) approach.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf

Implements the core of the nucleic acid folding algorithm. It calculates the
thermodynamically most stable secondary structure for a given single-stranded
DNA or RNA sequence. The result is a vector of  [`Structure`](@ref) objects, each element representing a
distinct secondary structure (e.g., hairpin loop, stacked pair, bulge, interior loop, multibranch loop)
that contributes to the overall folded structure.

An optimization is applied where "isolated" base pairs (those not adjacent to other base pairs)
are penalized with a high energy cost (1600.0 kcal/mol) to speed up computation. A base pair (i,j)
is considered isolated if neither the pair (i-1, j+1) nor the pair (i+1, j-1) are complementary
according to the sequence's complementarity rules. This optimization is applied regardless of
sequence length.

# Arguments
- `seq::AbstractString`: DNA/RNA sequence
- `temp::Real`: The temperature (°C) at which the folding is performed (default: `37.0`).


# Examples
```jldoctest
julia> fold("CCAACCGGTTGG")
4-element Vector{SeqFold.Structure}:
    1   12   -1.8  STACK:CC/GG    
    2   11   -1.5  STACK:CA/GT    
    3   10   -1.0  STACK:AA/TT    
    4    9    3.5  HAIRPIN:AC/TG  

julia> fold("CCAACCGGTTGG", temp=70)
2-element Vector{SeqFold.Structure}:
    5   12   -1.6  STACK_DE:CC/GG 
    6   11    3.2  HAIRPIN:CG/GT  
```

# See also
[`dg`](@ref), [`dg_cache`](@ref), [`dot_bracket`](@ref), [`Structure`](@ref),
[`SeqFold.DNA_ENERGIES`](@ref), [`SeqFold.RNA_ENERGIES`](@ref)
"""
function fold(seq::AbstractString; temp::Real = 37.0)::Vector{Structure}    
    v_cache, w_cache = _cache(seq, temp)
    n = length(seq)
    return _traceback(1, n, v_cache, w_cache)
end

"""
    dg(seq; temp = 37.0) -> Float64
    dg(structures) -> Float64

Compute the minimum free energy (ΔG, kcal/mol⁻¹) of the secondary structure
predicted for a single-stranded nucleic acid sequence at a specified temperature.

The function is a thin wrapper around the more general `fold` routine, which
generates all energetically feasible secondary structures for the sequence.  
Only the sum of the free-energy contributions of the returned structures is
reported, rounded to two decimal places.

# Arguments
- `seq::AbstractString` – the nucleotide sequence to be folded;
- `temp::Real` – the temperature (°C) at which to perform the folding
  (default: `37.0`);
- `structures::Vector{Structure}` – result of [`fold`](@ref) function.

# Returns
- `ΔG::Float64` – the total free energy of the predicted structure,
  rounded to two decimal places.

# Examples
```jldoctest
julia> seq = "GCGCGCGCGCG";

julia> dg(seq)
-3.0

julia> dg(seq, temp=20)
-4.9

julia> structures = fold(seq);

julia> dg(structures)
-3.0
```

# See also
[`fold`](@ref), [`dg_cache`](@ref), [`SeqFold.DG_ENERGIES`](@ref)
"""
function dg(seq::AbstractString; temp::Real = 37.0)::Float64
    structs = fold(seq, temp=temp)
    ΔG = sum(s.e for s in structs)
    return round(ΔG, digits=2)
end

function dg(structures::Vector{SeqFold.Structure})::Float64
    ΔG = sum(s.e for s in structures)
    return round(ΔG, digits=2)
end

"""Fold a nucleic acid sequence and return the estimated ΔG of each (i,j) pairing.

Args:
    seq: The nucleic acid sequence to fold

Keyword Args:
    temp: The temperature to fold at

Returns:
    Cache: A Vector{Vector} where each [i][j] pairing corresponds to the
        minimum free energy between i and j
"""
function dg_cache(seq::AbstractString, temp::Real = 37.0)::Cache
    _, w_cache = _cache(seq, temp)
    cache::Cache = [[s.e for s in row] for row in w_cache]
    return cache
end

"""
    dot_bracket(seq, structs) -> String

Generate the dot-bracket notation representation of a predicted nucleic acid secondary structure.

# Arguments
- `seq::AbstractString`: The original nucleotide sequence that was folded.
- `structs::Vector{Structure}`: A vector of `Structure` objects describing the folded structure,
  typically obtained from the [`fold`](@ref) function.

# Examples
```jldoctest
julia> s = "AATTACGTTAC";

julia> dot_bracket(s, fold(s))
"((.....)).."

julia> seq2 = "GGGAGGTCAGCAAACCTGAACCTGTTGAGATGTTGACGTCAGGAAACCCT";

julia> structs2 = fold(seq2)
12-element Vector{SeqFold.Structure}:
    3   40   -1.3  STACK:GA/CT    
    4   39   -1.4  STACK:AGG/TGC  
    6   37   -1.5  STACK:GT/CA    
    7   36   -1.3  STACK:TC/AG    
    8   35   -1.5  STACK:CA/GT    
    9   34   -2.0  STACK:AGC/TTG  
   11   32   -1.5  STACK:CA/GT    
   12   31    2.8  INTERIOR_LOOP:3/1
   16   29   -1.3  STACK:CT/GA    
   17   28   -0.1  STACK:TGA/AGT  
   19   26   -1.0  STACK:AA/TT    
   20   25    3.5  HAIRPIN:AC/TG  

julia> dbn = dot_bracket(seq2, structs2)
"..((.((((.((...((.((....)).)).)).)))).)).........."

julia> length(dbn) == length(seq2)
true

julia> count(==('('), dbn)  # Count base pairs
12
```

# Notes
- The function only considers the base pairs explicitly listed in the `.ij` field of each `Structure`.
- It does not validate that the input structures are consistent or represent a physically possible configuration.
- Pseudoknots (non-nested base pairs) are not handled by this simple representation and will not be correctly shown
  if present in the input structures.

# See also
[`fold`](@ref), [`Structure`](@ref)
"""
function dot_bracket(seq::AbstractString, structs::Vector{Structure})
    n = length(seq)
    result = fill('.', n)
    for s in structs
        if length(s.ij) == 1
            i, j = only(s.ij)
            result[i] = '('
            result[j] = ')'
        end
    end
    return String(result)
end

"""Create caches for the w_cache and v_cache

The Structs is useful for gathering many possible energies
between a series of `(i,j)` combinations.

Args:
    seq: The sequence to fold

Keyword Args:
    temp: The temperature to fold at

Returns:
    `(w_cache, v_cache)` tuple
"""
function _cache(seq, temp)
    if !(-50 ≤ temp ≤ 150)
        throw(ArgumentError("Temperature in °C is outside reasonable range"))
    end
    temp_K = temp + 273.15

    seq, n, _, is_dna = _parse_input(seq, is_dna=nothing)

    emap = is_dna ? DNA_ENERGIES : RNA_ENERGIES

    v_cache = [fill(STRUCT_DEFAULT, n) for _ in 1:n]
    w_cache = [fill(STRUCT_DEFAULT, n) for _ in 1:n]

    _w!(seq, 1, n, temp_K, v_cache, w_cache, emap)

    return v_cache, w_cache
end

"""
Find and return the lowest free energy structure in Sij subsequence (1-based indices).
This function calculates and stores results in `w_cache` to avoid recomputation.
Figure 2B in Zuker and Stiegler, 1981.

Args:
    seq: The sequence being folded.
    i: The start index.
    j: The end index (1-based, inclusive).
    temp: The temperature in Kelvin.
    v_cache: Free energy cache for if i and j base pair. Stores `Structure`.
    w_cache: Free energy cache for lowest energy structure from i to j. Stores `Structure`.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Structure: The minimum free energy structure for the subsequence from i to j.
"""
function _w!(seq, i, j, temp, v_cache, w_cache, emap)::Structure

    
    if w_cache[i][j] != STRUCT_DEFAULT
        return w_cache[i][j]
    end

    if j - i < 4
        w_cache[i][j] = STRUCT_NULL
        return w_cache[i][j]
    end


    w1 = _w!(seq, i + 1, j, temp, v_cache, w_cache, emap)
    w2 = _w!(seq, i, j - 1, temp, v_cache, w_cache, emap)
    w3 = _v!(seq, i, j, temp, v_cache, w_cache, emap)

    w4 = STRUCT_NULL

    for k in (i + 1):(j - 2) 
        w4_test = _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, false)
        if Bool(w4_test) && w4_test.e < w4.e
            w4 = w4_test
        end
    end
    
    return w_cache[i][j] = _min_struct(w1, w2, w3, w4)
end

"""
Find, store and return the minimum free energy of the structure between i and j.
If i and j don't base pair, store and return INF.
See: Figure 2B of Zuker, 1981

Args:
    seq: The sequence being folded (uppercase).
    i: The start index.
    j: The end index (1-based, inclusive).
    temp: The temperature in Kelvin.
    v_cache: Free energy cache for if i and j base pair. Stores `Structure`.
    w_cache: Free energy cache for lowest energy structure from i to j. Stores `Structure`.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Structure: The minimum free energy structure possible between i and j.
"""
function _v!(seq, i, j, temp, v_cache, w_cache, emap)::Structure
    # --- Base Case: Check Cache ---
    if v_cache[i][j] != STRUCT_DEFAULT
        return v_cache[i][j]
    end

    # --- Base Case: Check Complementarity ---
    # The ends must basepair for V(i,j)
    if emap.COMPLEMENT[seq[i]] != seq[j]
        v_cache[i][j] = STRUCT_NULL
        return v_cache[i][j]
    end

    n = length(seq)

    # if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
    # heuristic for speeding this up
    # from https://www.ncbi.nlm.nih.gov/pubmed/10329189
    isolated_outer = true
    if i > 1 && j < n
        isolated_outer = emap.COMPLEMENT[seq[i - 1]] != seq[j + 1]
    end
    isolated_inner = emap.COMPLEMENT[seq[i + 1]] != seq[j - 1]

    if isolated_outer && isolated_inner
        v_cache[i][j] = Structure(1600.0) # kcal/mol
        return v_cache[i][j]
    end

    # --- Energy E1: Hairpin Loop ---
    pair_hairpin = _pair(seq, i, i + 1, j, j - 1)
    e1_energy = _hairpin(seq, i, j, temp, emap)
    e1 = Structure(e1_energy, "HAIRPIN:" * pair_hairpin)
    
    # Special case: Small hairpin (4 bases in loop, so j-i=5)
    if j - i == 4
        v_cache[i][j] = e1
        w_cache[i][j] = e1
        return v_cache[i][j]
    end

    # --- Energy E2: Stacking/Bulge/Interior Loop ---
    # Stacking region or bulge or interior loop; Figure 2A(2)
    # j-i=d>4; various pairs i',j' for j'-i'<d
    e2 = Structure(Inf)
    for i1 in (i + 1):(j - 3)
        for j1 in (i1 + 4):(j - 1)
            # i1 and j1 must match (base pair)
            if emap.COMPLEMENT[seq[i1]] != seq[j1]
                continue
            end

            pair = _pair(seq, i, i1, j, j1)
            pair_left = _pair(seq, i, i + 1, j, j - 1)
            pair_right = _pair(seq, i1 - 1, i1, j1 + 1, j1)
            pair_inner = (haskey(emap.NN, pair_left) || haskey(emap.NN, pair_right))

            stack = (i1 == i + 1) && (j1 == j - 1)
            bulge_left = i1 > i + 1
            bulge_right = j1 < j - 1

            e2_test_energy = Inf
            e2_test_type = ""

            if stack
                # it's a neighboring/stacking pair in a helix
                e2_test_energy = _stack(seq, i, i1, j, j1, temp, emap)
                e2_test_type = "STACK:" * pair
                if (i > 1 && j == n) || (i == 1 && j < n)
                     # there's a dangling end
                     e2_test_type = "STACK_DE:" * pair
                end
            elseif bulge_left && bulge_right && !pair_inner
                # it's an interior loop
                loop_size_left = i1 - i
                loop_size_right = j - j1
                e2_test_energy = _internal_loop(seq, i, i1, j, j1, temp, emap)
                e2_test_type = string("INTERIOR_LOOP:", loop_size_left,'/', loop_size_right)
                if loop_size_left == 2 && loop_size_right == 2
                    loop_left = seq[i:i1]
                    loop_right = seq[j1:j]
                    # technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = string("STACK:", loop_left, '/', reverse(loop_right))
                end
            elseif bulge_left && !bulge_right
                # it's a bulge on the left side
                e2_test_energy = _bulge(seq, i, i1, j, j1, temp, emap)
                loop_size_left = i1 - i
                e2_test_type = string("BULGE:", loop_size_left)
            elseif !bulge_left && bulge_right
                # it's a bulge on the right side
                e2_test_energy = _bulge(seq, i, i1, j, j1, temp, emap)
                loop_size_right = j - j1
                e2_test_type = string("BULGE:", loop_size_right)
            else
                # it's basically a hairpin, only outside bp match, or other non-standard case
                continue
            end

            e2_test_energy += _v!(seq, i1, j1, temp, v_cache, w_cache, emap).e
            if e2_test_energy > -Inf && e2_test_energy < e2.e
                e2 = Structure(e2_test_energy, e2_test_type, [(i1, j1)])
            end
        end
    end

    # --- Energy E3: Multibranch Loop ---
    e3 = STRUCT_NULL
    if !isolated_outer || i == 1 || j == n
        for k in (i + 1):(j - 2)
            e3_test = _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, true)
            if Bool(e3_test) && e3_test.e < e3.e
                e3 = e3_test
            end
        end
    end

    return v_cache[i][j] = _min_struct(e1, e2, e3)
end

"""Return a stack representation, a key for the NN maps

Args:
    s: Sequence being folded
    i: leftmost index
    i1: index to right of i
    j: rightmost index
    j1: index to left of j

Returns:
    string representation of the pair
"""
function _pair(s, i, i1, j, j1)::String
    g(i) = get(s, i, '.')
    return string(g(i), g(i1), '/', g(j), g(j1))
end

"""Return the struct with the lowest free energy that isn't -inf (undef)

Args:
    structs: Structures being compared

Returns:
    struct: The min free energy structure
"""
function _min_struct(structs...)
    str = STRUCT_NULL
    for s in structs
        if s.e > -Inf && s.e < str.e
            str = s
        end
    end
    return str
end

"""Find the free energy given delta h, s and temp

Args:
    `d_h`: The enthalpy increment in kcal / mol
    `d_s`: The entropy increment in cal / mol
    `temp_K`: The temperature in Kelvin

Returns:
    The free energy increment in kcal / (mol x K)
"""
function _d_g(d_h, d_s, temp_K)
    d_h - temp_K * d_s / 1000.0
end

"""Estimate the free energy of length `query_len` based on one of length `known_len`.

The Jacobson-Stockmayer entry extrapolation formula is used
for bulges, hairpins, etc that fall outside the 30nt upper limit
for pre-calculated free-energies. See SantaLucia and Hicks (2004).

Args:
    query_len: Length of element without known free energy value
    known_len: Length of element with known free energy value (d_g_x)
    d_g_x: The free energy of the element known_len
    temp_K: Temperature in Kelvin

Returns:
    Free energy for a structure of length `query_len`
"""
function _j_s(query_len, known_len, d_g_x, temp_K)::Float64
    gas_constant = 1.9872e-3  # kcal/mol·K
    return d_g_x + 2.44 * gas_constant * temp_K * log(query_len / known_len)
end

"""
Get the free energy for a stack.
Using the indexes i and j, check whether it's at the end of
the sequence or internal. Then check whether it's a match
or mismatch, and return.
Two edge-cases are terminal mismatches and dangling ends.
The energy of a dangling end is added to the energy of a pair
where i XOR j is at the sequence's end.

Args:
    seq: The full folding sequence (uppercase).
    i: The start index on left side of the pair/stack.
    i1: The index to the right of i.
    j: The end index on right side of the pair/stack.
    j1: The index to the left of j.
    temp: Temperature in Kelvin.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Float64: The free energy of the NN pairing.
"""
function _stack(seq, i, i1, j, j1, temp, emap)::Float64
    n = length(seq)
    if any(>(n), (i, i1, j, j1))
        return 0.0
    end

    pair = _pair(seq, i, i1, j, j1)
    if any(<(1), (i, i1, j, j1)) # it's a dangling end
        d_h, d_s = emap.DE[pair]
        return _d_g(d_h, d_s, temp)
    end
    
    if i > 1 && j < n # it's internal
        d_h, d_s = haskey(emap.NN, pair) ? emap.NN[pair] : emap.INTERNAL_MM[pair]
        return _d_g(d_h, d_s, temp)
    end
    
    if i == 1 && j == n # it's terminal
        d_h, d_s = haskey(emap.NN, pair) ? emap.NN[pair] : emap.TERMINAL_MM[pair]
        return _d_g(d_h, d_s, temp)
    end

    if i > 1 && j == n # it's dangling on left
        d_h, d_s = haskey(emap.NN, pair) ? emap.NN[pair] : emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = string(seq[i - 1], seq[i], "/.", seq[j])
        if haskey(emap.DE, pair_de)
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        end
        return d_g
    end
    
    if i == 1 && j < n # it's dangling on right
        d_h, d_s = haskey(emap.NN, pair) ? emap.NN[pair] : emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = string('.', seq[i], '/', seq[j + 1], seq[j])
        if haskey(emap.DE, pair_de)
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        end
        return d_g
    end
    return 0.0
end

"""
Calculate the free energy of a hairpin.

Args:
    seq: The sequence we're folding (uppercase).
    i: The index of the start of hairpin.
    j: The index of the end of hairpin (1-based, inclusive).
    temp: Temperature in Kelvin.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Float64: The free energy increment from the hairpin structure.
"""
function _hairpin(seq, i, j, temp, emap)::Float64
    # --- Length Check ---
    if j - i < 4
        return Inf
    end

    hairpin = SubString(seq, i, j)
    hairpin_len = length(hairpin) - 2
    pair = _pair(seq, i, i + 1, j, j - 1)

    # --- Terminal Pair Check ---
    if emap.COMPLEMENT[first(hairpin)] != last(hairpin)
        # not known terminal pair, nothing to close "hairpin"
        throw(RuntimeError("Hairpin terminal pair does not match complement rules."))
    end

    # --- Initialize Energy ---
    d_g::Float64 = 0.0

    if !isnothing(emap.TRI_TETRA_LOOPS) && haskey(emap.TRI_TETRA_LOOPS, hairpin)
        # it's a pre-known hairpin with known value
        d_h, d_s = emap.TRI_TETRA_LOOPS[hairpin]
        d_g = _d_g(d_h, d_s, temp)
    end

    # --- Size-Based Penalty ---
    if 1 <= hairpin_len <= 30
        d_h, d_s = emap.HAIRPIN_LOOPS[hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else
        # Access the 30th element of the tuple for extrapolation
        d_h, d_s = emap.HAIRPIN_LOOPS[30]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, 30, d_g_inc, temp)
    end

    # Add penalty if loop is long enough (> 3 unpaired bases) and pair exists
    if hairpin_len > 3 && haskey(emap.TERMINAL_MM, pair)
         d_h, d_s = emap.TERMINAL_MM[pair]
         d_g += _d_g(d_h, d_s, temp)
    end

    # --- Special AT Closing Penalty (SantaLucia 2004) ---
    # Add penalty if loop length is exactly 3 and closing pair involves 'A'
    if hairpin_len == 3 && (first(hairpin) == 'A' || last(hairpin) == 'A')
        d_g += 0.5 # kcal/mol
    end

    return d_g
end

"""
Calculate the free energy associated with a bulge.

Args:
    seq: The full folding DNA/RNA sequence (uppercase).
    i: The start index of the bulge.
    i1: The index to the right of i.
    j: The end index of the bulge.
    j1: The index to the left of j.
    temp: Temperature in Kelvin.
    emap: Map to DNA/RNA energies (`Energies` struct).

Returns:
    Float64: The increment in free energy from the bulge.
"""
function _bulge(seq, i, i1, j, j1, temp, emap)::Float64
    # --- Calculate Bulge Loop Length ---
    loop_len = max(i1 - i - 1, j - j1 - 1)
    if loop_len <= 0
       throw(RuntimeError("Bulge loop length must be positive."))
    end

    # --- Initialize Energy ---
    d_g = 0.0

    # add penalty based on size
    if 1 <= loop_len <= 30
        d_h, d_s = emap.BULGE_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else
        # it's too large for pre-calculated list, extrapolate
        d_h, d_s = emap.BULGE_LOOPS[30]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g, temp)
    end

    if loop_len == 1
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        pair = _pair(seq, i, i1, j, j1)
        @assert pair in keys(emap.NN)
        d_g += _stack(seq, i, i1, j, j1, temp, emap)
    end
    # penalize AT terminal bonds
    if any(k->(seq[k]=='A'), (i, i1, j, j1))
        d_g += 0.5
    end
    return d_g
end

"""
Calculate the free energy of an internal loop.

The first and last bp of both left and right sequences
are not themselves parts of the loop, but are the terminal
bp on either side of it. They are needed for when there's
a single internal looping bp (where just the mismatching
free energies are used).

Note that both left and right sequences are in 5' to 3' direction.

This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004.

Args:
    seq: The sequence we're folding (uppercase).
    i: The index of the start of structure on left side.
    i1: The index to the right of i.
    j: The index of the end of structure on right side.
    j1: The index to the left of j.
    temp: Temperature in Kelvin.
    emap: Dictionary mapping to energies for DNA/RNA (`Energies` struct).

Returns:
    Float64: The free energy associated with the internal loop.
"""
function _internal_loop(seq, i, i1, j, j1, temp, emap)::Float64
    # --- Calculate Loop Sizes ---
    loop_left = i1 - i - 1
    loop_right = j - j1 - 1
    loop_len = loop_left + loop_right

    # --- Validate Loop Sizes ---
    if loop_left < 1 || loop_right < 1
        throw(RuntimeError("Internal loop sides must each have at least one unpaired base."))
    end

    # --- Special Case: Single Base Mismatch (1x1 Internal Loop) ---
    if loop_left == 1 && loop_right == 1
        # single bp mismatch, sum up the two single mismatch pairs
        mm_left = _stack(seq, i, i1, j, j1, temp, emap)
        mm_right = _stack(seq, i1 - 1, i1, j1 + 1, j1, temp, emap)
        return mm_left + mm_right
    end

    # --- Initialize Energy ---
    d_g::Float64 = 0.0

    if 1 <= loop_len <= 30
        d_h, d_s = emap.INTERNAL_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else
        # Access the 30th element of the tuple for extrapolation
        d_h, d_s = emap.INTERNAL_LOOPS[30]
        d_g_base = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g_base, temp)
    end

    # --- Asymmetry Penalty ---
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry # kcal/mol penalty per asymmetric base

    # --- Terminal Mismatch Penalties ---
    # Apply penalty based on the mismatching pairs on either side of the loop

    # Left side
    pair_left_mm = _pair(seq, i, i + 1, j, j - 1)
    d_h, d_s = emap.TERMINAL_MM[pair_left_mm]
    d_g += _d_g(d_h, d_s, temp)

    # Right side
    pair_right_mm = _pair(seq, i1 - 1, i1, j1 + 1, j1)
    d_h, d_s = emap.TERMINAL_MM[pair_right_mm]
    d_g += _d_g(d_h, d_s, temp)

    return d_g
end

"""
Calculate a multi-branch energy penalty using a linear formula.

From Jaeger, Turner, and Zuker, 1989.
Found to be better than logarithmic in Ward, et al. 2017.

Args:
    seq: The sequence being folded (uppercase).
    i: The left starting index.
    k: The mid-point in the search.
    j: The right ending index.
    temp: Folding temp in Kelvin.
    v_cache: Cache of energies where V(i,j) bond. Stores `Structure`.
    w_cache: Cache of min energy of substructures between W(i,j). Stores `Structure`.
    emap: Map to DNA/RNA energies (`Energies` struct).
    helix: Whether this multibranch is enclosed by a helix (V(i,j) forms a base pair).

Returns:
    Structure: A multi-branch structure.
"""
function _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, helix=false)::Structure
    left = helix ?
        _w!(seq, i + 1, k, temp, v_cache, w_cache, emap) :
        _w!(seq, i, k, temp, v_cache, w_cache, emap)
    right = helix ?
        _w!(seq, k + 1, j - 1, temp, v_cache, w_cache, emap) :
        _w!(seq, k + 1, j, temp, v_cache, w_cache, emap)

    if !Bool(left) || !Bool(right)
        return STRUCT_NULL
    end

    branches = Tuple{Int, Int}[]

    function __add_branch!(s::Structure)
        if !Bool(s) || isempty(s.ij)
            return
        end
        if length(s.ij) == 1
            push!(branches, s.ij[1])
            return
        end
        for (i1, j1) in s.ij
            sub_struct = _w!(seq, i1, j1, temp, v_cache, w_cache, emap)
            __add_branch!(sub_struct)
        end
    end

    __add_branch!(left)
    __add_branch!(right)

    # --- Validate Number of Branches ---
    # A multibranch must have at least 2 distinct branches.
    if length(branches) < 2
        return STRUCT_NULL
    end

    # --- Include Enclosing Helix (if applicable) ---
    # If this multibranch is enclosed by a helix (i,j), add it to the list.
    if helix
        push!(branches, (i, j))
    end

    # --- Calculate Energy Components ---
    branches_count = length(branches)
    unpaired = 0
    e_sum = 0.0

    for (index, (i2, j2)) in enumerate(branches)
        _, j1 = branches[mod1(index - 1, branches_count)]
        i3, j3 = branches[mod1(index + 1, branches_count)]

        # --- Dangling End Energy Calculation ---
        unpaired_left = 0
        unpaired_right = 0
        de = 0.0

        # Special handling for the last branch in the non-helix case
        if index == length(branches) && !helix
            # Do nothing for the last branch if not enclosed by a helix
        elseif (i3, j3) == (i, j)
            unpaired_left = i2 - j1 - 1
            unpaired_right = j3 - j2 - 1

            if unpaired_left != 0 && unpaired_right != 0
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elseif unpaired_right != 0
                de = _stack(seq, -1, i2, j2 + 1, j2, temp, emap)
                if unpaired_right == 1
                    de = min(_stack(seq, i3, -1, j3, j3 - 1, temp, emap), de)
                end
            end
        elseif (i2, j2) == (i, j)
            unpaired_left = j2 - j1 - 1
            unpaired_right = i3 - i2 - 1

            if unpaired_left != 0 && unpaired_right != 0
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elseif unpaired_right != 0
                de = _stack(seq, i2, i2 + 1, j2, -1, temp, emap)
                if unpaired_right == 1
                    de = min(_stack(seq, i3 - 1, i3, -1, j3, temp, emap), de)
                end
            end
        else
            unpaired_left = i2 - j1 - 1
            unpaired_right = i3 - j2 - 1

            if unpaired_left != 0 && unpaired_right != 0
                de = _stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap)
            elseif unpaired_right != 0
                de = _stack(seq, -1, i2, j2 + 1, j2, temp, emap)
                if unpaired_right == 1
                    de = min(_stack(seq, i2 - 1, i2, j2 + 1, j2, temp, emap), de)
                end
            end
        end

        e_sum += de
        unpaired += unpaired_right
        @assert unpaired_right >= 0

        if (i2, j2) != (i, j)
            e_sum += _w!(seq, i2, j2, temp, v_cache, w_cache, emap).e
        end
    end
    @assert unpaired >= 0

    a, b, c, d = emap.MULTIBRANCH

    # Calculate the multibranch loop initiation/penalty energy
    e_multibranch = unpaired == 0 ? a + d : a + b * length(branches) + c * unpaired
    e = e_multibranch + e_sum

    if helix
        pop!(branches)
    end

    return Structure(e, "BIFURCATION:" * string(unpaired) * "n/" * string(branches_count) * "h", branches)
end

"""
Traceback thru the V(i,j) and W(i,j) caches to find the structure.

For each step, get to the lowest energy W(i,j) within that block.
Store the structure in W(i,j).
Inc i and j.
If the next structure is viable according to V(i,j), store as well.
Repeat.

Args:
    i: The leftmost index to start searching in.
    j: The rightmost index to start searching in.
    v_cache: Cache of energies where i and j bond. Stores `Structure`.
    w_cache: Cache of energies/sub-structures between or with i and j. Stores `Structure`.

Returns:
    Vector{Structure}: A list of Structs in the final secondary structure.
"""
function _traceback( i, j, v_cache, w_cache)::Vector{Structure}
    s_w = w_cache[i][j]
    if !occursin("HAIRPIN", s_w.desc)
        while w_cache[i + 1][j] == s_w
            i += 1
        end
        while w_cache[i][j - 1] == s_w
            j -= 1
        end
    end

    structs::Vector{Structure} = Structure[]
    while true
        s_v = v_cache[i][j]

        # --- Handle Multibranch Structures ---
        if length(s_w.ij) > 1
            s_v = s_w
        end

        push!(structs, with_ij(s_v, [(i, j)]))

        # --- Case 1: Multibranch ---
        if length(s_v.ij) > 1
            e_sum = 0.0
            structs = _trackback_energy(structs)
            branches = Structure[]
            for (i1, j1) in s_v.ij
                tb = _traceback(i1, j1, v_cache, w_cache)
                if !isempty(tb) && !isempty(tb[1].ij)
                    i2, j2 = tb[1].ij[1]
                    e_sum += w_cache[i2][j2].e
                    append!(branches, tb)
                end
            end
            last_struct = last(structs)
            corrected_energy = round(last_struct.e - e_sum, digits=1)
            structs[end] = Structure(corrected_energy, last_struct.desc, collect(last_struct.ij))
            return vcat(structs, branches)
        end

        # --- Case 2: Single Structure (Stack, Bulge, etc.) ---
        if length(s_v.ij) == 1
            i, j = only(s_v.ij)
            continue
        end

        # --- Case 3: Hairpin (End of Structure) ---
        return _trackback_energy(structs)
    end
end

"""
Add energy to each structure, based on how it's W(i,j) differs from the one after.

Args:
    structs: The structures for whom energy is being calculated.

Returns:
    Vector{Structure}: Structures in the folded DNA with energy.
"""
function _trackback_energy(structs)::Vector{Structure}
    n = length(structs)
    structs_e = Vector{Structure}(undef, n)

    @inbounds for (i, str) in enumerate(structs)
        e_next = i == n ? 0.0 : structs[i+1].e
        e_corrected = round(str.e - e_next, digits=1)
        corrected_struct = Structure(e_corrected, str.desc, copy(str.ij))
        structs_e[i] = corrected_struct
    end

    return structs_e
end