    """Fold the DNA sequence and return the lowest free energy score.

Based on the approach described in:
Zuker and Stiegler, 1981
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf

If the sequence is 50 or more bp long, "isolated" matching bp
are ignored in V(i,j). This is based on an approach described in:
Mathews, Sabina, Zuker and Turner, 1999
https://www.ncbi.nlm.nih.gov/pubmed/10329189

Args:
    seq: The sequence to fold

Keyword Args:
    temp: The temperature the fold takes place in, in Celcius

Returns:
    A vector of structures. Stacks, bulges, hairpins, etc.
"""
function fold(seq::AbstractString; temp::Real = 37.0)::Vector{Struct}    
    v_cache, w_cache = _cache(seq, temp)
    n = length(seq)
    return _traceback(1, n, v_cache, w_cache)
end

"""Fold the sequence and return just the delta G of the structure

Args:
    seq: The sequence to fold

Keyword Args:
    temp: The temperature to fold at

Returns:
    Minimum free energy of the folded sequence
"""
function dg(seq::AbstractString; temp::Real = 37.0)::Float64
    structs = fold(seq, temp=temp)
    ΔG = sum(s.e for s in structs)
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

"""Get the dot bracket notation for a secondary structure.

Args:
    structs: A list of structs, usually from the fold function

Returns:
    Dot bracket notation of the secondary structure
"""
function dot_bracket(seq::AbstractString, structs::Vector{Struct})
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
    end

    if all(b in "AUCG" for b in unique_bases)
        is_dna = false
    elseif any(b ∉ "ATGC" for b in unique_bases)
        unknown_bases = filter(b -> b ∉ "ATUGC", unique_bases)
        throw(ArgumentError("Unknown base(s): $unknown_bases. Only DNA/RNA foldable."))
    end

    emap = is_dna ? DNA_ENERGIES : RNA_ENERGIES

    n = length(seq_str)
    v_cache = [fill(STRUCT_DEFAULT, n) for _ in 1:n]
    w_cache = [fill(STRUCT_DEFAULT, n) for _ in 1:n]

    _w!(seq_str, 1, n, temp_K, v_cache, w_cache, emap)

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
    v_cache: Free energy cache for if i and j base pair. Stores `Struct`.
    w_cache: Free energy cache for lowest energy structure from i to j. Stores `Struct`.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Struct: The minimum free energy structure for the subsequence from i to j.
"""
function _w!(seq, i, j, temp, v_cache, w_cache, emap)::Struct

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _w! ($(i-1), $(j-1))")

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
    v_cache: Free energy cache for if i and j base pair. Stores `Struct`.
    w_cache: Free energy cache for lowest energy structure from i to j. Stores `Struct`.
    emap: Energy map for DNA/RNA (`Energies` struct).

Returns:
    Struct: The minimum free energy structure possible between i and j.
"""
function _v!(seq, i, j, temp, v_cache, w_cache, emap)::Struct

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _v! ($(i-1), $(j-1))")


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
        v_cache[i][j] = Struct(1600.0) # kcal/mol
        return v_cache[i][j]
    end

    # --- Energy E1: Hairpin Loop ---
    pair_hairpin = _pair(seq, i, i + 1, j, j - 1)
    e1_energy = _hairpin(seq, i, j, temp, emap)
    e1 = Struct(e1_energy, "HAIRPIN:" * pair_hairpin)
    
    # Special case: Small hairpin (4 bases in loop, so j-i=5)
    if j - i == 4
        v_cache[i][j] = e1
        w_cache[i][j] = e1
        return v_cache[i][j]
    end

    # --- Energy E2: Stacking/Bulge/Interior Loop ---
    # Stacking region or bulge or interior loop; Figure 2A(2)
    # j-i=d>4; various pairs i',j' for j'-i'<d
    e2 = Struct(Inf)
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
                loop_size_left = i1 - i - 1
                loop_size_right = j - j1 - 1
                e2_test_energy = _internal_loop(seq, i, i1, j, j1, temp, emap)
                e2_test_type = string("INTERIOR_LOOP:", loop_size_left,'/', loop_size_right)
                if loop_size_left == 1 && loop_size_right == 1
                    loop_left = seq[i:i1]
                    loop_right = seq[j1:j]
                    # technically an interior loop of 1. really 1bp mismatch
                    e2_test_type = string("STACK:", loop_left, '/', reverse(loop_right))
                end
            elseif bulge_left && !bulge_right
                # it's a bulge on the left side
                e2_test_energy = _bulge(seq, i, i1, j, j1, temp, emap)
                loop_size_left = i1 - i - 1
                e2_test_type = string("BULGE:", loop_size_left)
            elseif !bulge_left && bulge_right
                # it's a bulge on the right side
                e2_test_energy = _bulge(seq, i, i1, j, j1, temp, emap)
                loop_size_right = j - j1 - 1
                e2_test_type = string("BULGE:", loop_size_right)
            else
                # it's basically a hairpin, only outside bp match, or other non-standard case
                continue
            end

            e2_test_energy += _v!(seq, i1, j1, temp, v_cache, w_cache, emap).e
            if e2_test_energy > -Inf && e2_test_energy < e2.e
                e2 = Struct(e2_test_energy, e2_test_type, [(i1, j1)])
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
    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _pair ($(i-1), $(i1-1), $(j-1), $(j1-1))")

    g(i) = get(s, i, '.')
    return string(
        g(i),
        g(i1),
        '/',
        g(j),
        g(j1)
)
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
    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _d_g ($d_h, $d_s)")

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
    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _j_s ($query_len, $known_len, $d_g_x)")


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
    
    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _stack ($(i-1), $(i1-1), $(j-1), $(j1-1))")


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

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _hairpin ($(i-1), $(j-1))")
    
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

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _bulge ($(i-1), $(i1-1), $(j-1), $(j1-1))")
    
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
    
    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _internal_loop ($(i-1), $(i1-1), $(j-1), $(j1-1))")

    n = length(seq)

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
    v_cache: Cache of energies where V(i,j) bond. Stores `Struct`.
    w_cache: Cache of min energy of substructures between W(i,j). Stores `Struct`.
    emap: Map to DNA/RNA energies (`Energies` struct).
    helix: Whether this multibranch is enclosed by a helix (V(i,j) forms a base pair).

Returns:
    Struct: A multi-branch structure.
"""
function _multi_branch(seq, i, k, j, temp, v_cache, w_cache, emap, helix=false)::Struct

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _multi_branch ($(i-1), $(k-1), $(j-1))")

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

    function __add_branch!(s::Struct)
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

    return Struct(e, "BIFURCATION:" * string(unpaired) * "n/" * string(branches_count) * "h", branches)
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
    v_cache: Cache of energies where i and j bond. Stores `Struct`.
    w_cache: Cache of energies/sub-structures between or with i and j. Stores `Struct`.

Returns:
    Vector{Struct}: A list of Structs in the final secondary structure.
"""
function _traceback( i, j, v_cache, w_cache)::Vector{Struct}

    # global COUNTER
    # COUNTER += 1
    # println("$COUNTER: _traceback ($(i-1), $(j-1))")

    s_w = w_cache[i][j]
    if !occursin("HAIRPIN", s_w.desc)
        while w_cache[i + 1][j] == s_w
            i += 1
        end
        while w_cache[i][j - 1] == s_w
            j -= 1
        end
    end

    structs::Vector{Struct} = Struct[]
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
            branches = Struct[]
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
            structs[end] = Struct(corrected_energy, last_struct.desc, collect(last_struct.ij))
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
    Vector{Struct}: Structures in the folded DNA with energy.
"""
function _trackback_energy(structs)::Vector{Struct}
    n = length(structs)
    structs_e = Vector{Struct}(undef, n)

    @inbounds for (i, str) in enumerate(structs)
        e_next = i == n ? 0.0 : structs[i+1].e
        e_corrected = round(str.e - e_next, digits=1)
        corrected_struct = Struct(e_corrected, str.desc, copy(str.ij))
        structs_e[i] = corrected_struct
    end

    return structs_e
end