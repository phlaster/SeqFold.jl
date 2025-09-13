function _parse_input(s1::AbstractString, s2::AbstractString)
    seq1, n1, gc_count1 = _parse_input(s1)
    seq2, n2, _ = _parse_input(s2)

    if n1 != n2
        throw(DimensionMismatch("Length mismatch between sequences: $n1 vs $n2"))
    end

    return seq1, seq2, n1, gc_count1
end

function _parse_input(seq::AbstractString; is_dna=true)
    n = length(seq)
    if n < 2
        throw(ArgumentError("Sequences lengths have to be at least 2bp"))
    end
    
    seq = uppercase(seq)

    if !isnothing(is_dna)
        gc_count = 0
        correct_bases = is_dna ? "ACGT" : "ACGU"
        for c in seq
            if !in(c, correct_bases)
                throw(ArgumentError("Incorrect character '$c' in $(is_dna ? "DNA" : "RNA") sequence"))
            end
            gc_count += Int(c == 'G' || c == 'C')
        end

        return seq, n, gc_count, is_dna
    end

    unique_bases = Set(seq)
    is_dna = true

    if 'U' in unique_bases && 'T' in unique_bases
        throw(ArgumentError("Both T and U found in sequence. Provide DNA or RNA, not both."))
    end

    if all(b in "AUCG" for b in unique_bases)
        is_dna = false
    elseif any(b ∉ "ATGC" for b in unique_bases)
        unknown_bases = filter(b -> b ∉ "ATUGC", unique_bases)
        throw(ArgumentError("Unknown base(s): $unknown_bases. Only DNA/RNA foldable."))
    end
    return seq, n, NaN, is_dna
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

"""
    SeqFold.complement(seq)

Compute the Watson-Crick complement of a DNA sequence.

This function converts each base in the input sequence to its complementary base,
handling both uppercase and lowercase inputs. Non-standard bases (anything other than `'A'`, `'T'`, `'C'`, `'G'`)
are converted to a period (`'.'`).

# Arguments
- `seq::AbstractString | AbstractVector{UInt8}`: A DNA sequence as either a string or raw byte representation.

# Returns
A new string (if input was a string) or vector of bytes (if input was bytes) representing the complementary DNA sequence.

# Examples
```jldoctest
julia> SeqFold.complement("ACGT")
"TGCA"

julia> SeqFold.complement("ACGTacgt")
"TGCATGCA"

julia> SeqFold.complement("ACGTN")
"TGCA."

julia> SeqFold.complement("XYZ")
"..."

julia> SeqFold.complement(b"ACGT")
4-element Vector{UInt8}:
 0x54
 0x47
 0x43
 0x41

julia> String(SeqFold.complement(b"ACGT"))
"TGCA"

julia> s = "ACACTAC"; SeqFold.complement(SeqFold.complement(s)) == s
true
```

# Notes
- The function normalizes all output to uppercase bases.
- For string inputs, the output is a string; for byte vector inputs, the output is a byte vector.
- This is specifically designed for DNA sequences (not RNA, which would use U instead of T).
- The complement of a complement returns the original sequence (except for non-standard bases which become '.').
"""
function complement(seq::AbstractString)
    out_bytes = complement(codeunits(seq))
    return String(out_bytes)
end