"""
    AbstractKmerIterator{A <: Alphabet, K}

Iterates `Kmer{A, K}`.
Functions to implement:
* `Base.iterate`

Optional functions:
* `source_type`
* `load_source`
"""
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

"""
    load_source(x::AbstractKmerIterator)::source_type(typeof(x))

Get the data source from the kmer iterator.
"""
function load_source end

"""
usable_source(x::AbstractKmerIterator)

Convert the source object into whatever is used by the iterator protocol
"""
function usable_source(x::AbstractKmerIterator)::Union{BioSequence, AbstractVector{UInt8}}
    loaded = load_source(x)
    return if loaded isa Union{BioSequence, AbstractVector{UInt8}}
        loaded
    elseif loaded isa Union{String, SubString{String}}
        codeunits(loaded)
    else
        error("Does not know how to load data from source")
    end
end

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
    elseif sT <: ByteSource && A isa Union{FourBit, AminoAcidAlphabet}
        Base.HasLength()
    else
        Base.SizeUnknown()
    end
end

"""Trait object which based on static dispatch determines how to recode from
the encoding of the source sequence to the encoding of the kmer"""
abstract type RecodingScheme end

"We can copy the encoding straight from the source to the kmer"
struct Copyable <: RecodingScheme end

"We can copy the encoding, then bitshift to create 4-bit encoding"
struct TwoToFour <: RecodingScheme end

"We skip all symbols in source that contain unmappable symbols"
struct Skipping <: RecodingScheme end

"We can use `BioSequences.ascii_encode`"
struct AsciiEncode <: RecodingScheme end

"The source is bytes, but we need our own encoding table,
since we must skip ambiguous nucleotides"
struct AsciiSkipping <: RecodingScheme end

"The source is a bytevector, and we have no static knowledge of efficient
conversion to the right encoding"
struct GenericBytes <: RecodingScheme end

"Generic fallback when the source is a `BioSequence`"
struct GenericAlphabet <: RecodingScheme end

function RecodingScheme(::Type{T})::RecodingScheme where {T <: AbstractKmerIterator}
    A = Alphabet(eltype(T))
    sT = source_type(T)
    return if sT <: BioSequence
        As = Alphabet(sT)
        if As == A
            Copyable()
        elseif As isa TwoBit && A isa TwoBit
            Copyable()
        elseif As isa FourBit && A isa FourBit
            Copyable()
        elseif As isa FourBit && A isa TwoBit
            Skipping()
        elseif As isa TwoBit && A isa FourBit
            TwoToFour()
        else
            GenericAlphabet()
        end
    elseif sT <: ByteSource
        codetype = BioSequences.codetype(A)
        return if codetype isa BioSequences.AsciiAlphabet
            if A isa TwoBit
                AsciiSkipping()
            else
                AsciiEncode()
            end
        else
            return GenericBytes()
        end
    else
        error("Cannot determine recoding scheme of iterator")
    end
end

@noinline throw_bad_byte_error(b::UInt8) =
    error("Cannot interpret byte $(repr(b)) as nucleotide")

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

const FOURBIT_COMPLEMENT_LUT = let
    v = fill(0x00, 16)
    for i in alphabet(DNA)
        v[reinterpret(UInt8, i) + 0x01] = reinterpret(UInt8, complement(i))
    end
    Tuple(v)
end

"Extract a full kmer at a given index of a sequence.
Note: These methods don't do any bounds checking"
function extract end
# TODO: Use extract elsewhere in this code base, e.g. kmer from string instantiation?

@inline function extract(
    ::GenericAlphabet,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    length(seq) < ksize(T) && return nothing
    data = zero_tuple(T)
    A = Alphabet(T)
    bps = BioSequences.bits_per_symbol(A)
    @inbounds for i in 1:ksize(T)
        symbol = seq[i]
        encoding = UInt(BioSequences.encode(A, convert(eltype(A), symbol)))::UInt
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end

@inline function extract(
    ::TwoToFour,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    length(seq) < ksize(T) && return nothing
    data = zero_tuple(T)
    for i in 1:ksize(T)
        encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, i)))
        (_, data) = leftshift_carry(data, 4, encoding)
    end
    T(unsafe, data)
end

@inline function extract(
    ::Copyable,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(seq))
    for i in from_index:(from_index + ksize(T) - 1)
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end

@inline function extract(
    ::AsciiEncode,
    ::Type{T},
    seq::AbstractVector{UInt8},
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(kT))
    @inbounds for i in from_index:(from_index + ksize(T) - 1)
        encoding = BioSequences.ascii_encode(Alphabet(T), seq[i])
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end

@inline function extract(
    ::GenericBytes,
    ::Type{T},
    seq::AbstractVector{UInt8},
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(T))
    @inbounds for i in 1:ksize(T)
        char = reinterpret(Char, (seq[i] % UInt32) << 24)
        symbol = eltype(T)(char)
        encoding = UInt(BioSequences.encode(Alphabet(T), symbol))::UInt
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end
