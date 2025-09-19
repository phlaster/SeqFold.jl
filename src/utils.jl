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
    UInt8('N')

const DNA_COMP_TABLE = collect(ntuple(i->_comp_dna(UInt8(i-1)), Val(256)))

@inline function _complement_bytes!(out::Vector{UInt8}, seq::AbstractVector{UInt8}, table)
    n = length(seq)
    @inbounds @simd for i in 1:n
        out[i] = table[ seq[i] + 0x01 ]
    end
    return out
end

function _complement_bytes(nucleotide::UInt8, table)
    table[nucleotide + 0x01]
end 

function _complement_bytes(seq::AbstractVector{UInt8}, table)
    out = Vector{UInt8}(undef, length(seq))
    return _complement_bytes!(out, seq, table)
end

function complement(nucleotide::AbstractChar; table=DNA_COMP_TABLE)
    Char(_complement_bytes(UInt8(nucleotide), table))
end


"""
    SeqFold.complement(::AbstractString) -> String
    SeqFold.complement(::AbstractChar) -> Char

Compute the Watson-Crick all uppercase complement of a DNA sequence.
Non-standard bases (anything other than `'A'`, `'T'`, `'C'`, `'G'`)
are converted to `'N'`.

# Examples
```jldoctest
julia> SeqFold.complement("ACGT")
"TGCA"

julia> SeqFold.complement("ACGTacgt")
"TGCATGCA"

julia> SeqFold.complement("ACGTN")
"TGCAN"

julia> SeqFold.complement("XYZ")
"NNN"

julia> SeqFold.complement('C')
'G': ASCII/Unicode U+0047 (category Lu: Letter, uppercase)
```

# See also
[`SeqFold.revcomp`](@ref)
"""
function complement(seq::AbstractString; table=DNA_COMP_TABLE)
    cu = codeunits(seq)
    if length(cu) != length(seq)
        throw(ErrorException("Some characters in $seq occupy >1 codeunits. ASCII characters only!"))
    end
    out_bytes = _complement_bytes(cu, table)
    return String(out_bytes)
end


@inline function _revcomp_bytes!(out::Vector{UInt8}, seq::AbstractVector{UInt8}, table)
    n = length(seq)
    @inbounds @simd for i in 1:n
        out[i] = table[ seq[n - i + 1] + 0x01 ]
    end
    return out
end

function revcomp(nucleotide::AbstractChar; table=DNA_COMP_TABLE)
    complement(nucleotide; table=table)
end


function _revcomp_bytes(seq::AbstractVector{UInt8}, table)
    out = Vector{UInt8}(undef, length(seq))
    return _revcomp_bytes!(out, seq, table)
end


"""
    SeqFold.revcomp(::AbstractString) -> String
    SeqFold.revcomp(::AbstractChar) -> Char

Compute the Watson-Crick all uppercase reverse complement of a DNA sequence.
Non-standard bases (anything other than `'A'`, `'T'`, `'C'`, `'G'`)
are converted to `'N'`.

# Examples
```jldoctest
julia> SeqFold.revcomp("ACGT")
"ACGT"

julia> SeqFold.revcomp("ACGTacgt")
"ACGTACGT"

julia> SeqFold.revcomp("ACGTN")
"NACGT"

julia> SeqFold.revcomp("XYZ")
"NNN"

julia> SeqFold.revcomp('C')
'G': ASCII/Unicode U+0047 (category Lu: Letter, uppercase)
```

# See also
[`SeqFold.complement`](@ref)
"""
function revcomp(seq::AbstractString; table=DNA_COMP_TABLE)
    cu = codeunits(seq)
    if length(cu) != length(seq)
        throw(ErrorException("Some characters in $seq occupy >1 codeunits. ASCII characters only!"))
    end
    out_bytes = _revcomp_bytes(cu, table)
    return String(out_bytes)
end


"""
    SeqFold.gc_content(::AbstractString) -> Float64

Compute the GC-content of a DNA(/RNA) sequence. Does not validate character correctness.

# Examples
```jldoctest
julia> SeqFold.gc_content("ACGT")
0.5

julia> SeqFold.gc_content("GGC")
1.0

julia> SeqFold.gc_content("ABCD")
0.25
```

# See also
[`SeqFold.gc_cache`](@ref)
"""
function gc_content(seq::AbstractString)
    counter = 0
    @simd for c in seq
        if c in "GCgc"
            counter += 1
        end
    end
    return counter / length(seq)
end