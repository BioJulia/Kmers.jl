abstract type AbstractKmerIterator{A <: Alphabet, K} end

function Base.eltype(::Type{<:AbstractKmerIterator{A, K}}) where {A, K}
    Kmer{A, K, n_coding_elements(Kmer{A, K})}
end

Base.IteratorSize(::Type{<:AbstractKmerIterator}) = Base.SizeUnknown()

const FourBit = Union{DNAAlphabet{4}, RNAAlphabet{4}}
const TwoBit = Union{DNAAlphabet{2}, RNAAlphabet{2}}

@noinline throw_bad_byte_error(b::UInt8) = error("Cannot interpret byte $(repr(b)) as nucleotide")

const BYTE_LUT = let
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