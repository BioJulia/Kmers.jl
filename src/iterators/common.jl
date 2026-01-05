"""
    AbstractKmerIterator{A <: Alphabet, K}

Abstract type for kmer iterators. The element type is `Kmer{A, K, N}`,
with the appropriately derived N.

Functions to implement:
* `Base.iterate`
* `Base.length` or `Base.IteratorSize` if not `HasLength`
"""
abstract type AbstractKmerIterator{A <: Alphabet, K} end

function Base.eltype(::Type{<:AbstractKmerIterator{A, K}}) where {A, K}
    return Kmer{A, K, n_coding_elements(Kmer{A, K})}
end

function used_source(::AsciiEncode, s::AbstractString)
    return is_ascii(typeof(s)) ? codeunits(s) : s
end
used_source(::RecodingScheme, s) = s

const ASCII_SKIPPING_LUT = let
    v = fill(0xff, 256)
    for (i, s) in [(0, "Aa"), (1, "cC"), (2, "gG"), (3, "TtUu")], c in s
        v[UInt8(c) + 1] = i
    end
    for c in "-MRSVWYHKDBN"
        v[UInt8(c) + 1] = 0xf0
        v[UInt8(lowercase(c)) + 1] = 0xf0
    end
    Tuple(v)
end
