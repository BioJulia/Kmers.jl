abstract type AbstractKmerIterator{A <: Alphabet, K} end

function Base.eltype(::Type{<:AbstractKmerIterator{A, K}}) where {A, K}
    Kmer{A, K, n_coding_elements(Kmer{A, K})}
end

const FourBit = Union{DNAAlphabet{4}, RNAAlphabet{4}}
const TwoBit = Union{DNAAlphabet{2}, RNAAlphabet{2}}
const ByteSource = Union{String, SubString{String}, AbstractVector{UInt8}}

"""
    source_type(::Type{<:AbstractKmerIterator})::Type

Get the type of the data source that kmers are extracted from
"""
function source_type end

# We can't compute the length for
# * Unknown alphabets, or unknown source types
# * A source type that can contain ambiguous nucleotides,
#   while the kmer type does not. In this case, it's standard practise to
#   skip these symbols.
function Base.IteratorSize(::Type{T}) where {T <: AbstractKmerIterator}
    kT = eltype(T)
    sT = source_type(T)
    A = Alphabet(kT)
    return if sT <: BioSequence && Alphabet(sT) == A
        Base.HasLength()
    elseif sT <: BioSequence && Alphabet(sT) isa TwoBit && A isa FourBit
        Base.HasLength()
    elseif sT <: ByteSource && isa Union{FourBit, AminoAcidAlphabet}
        Base.HasLength()
    else
        Base.SizeUnknown()
    end
end

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