################################################
# Trait dispatch
################################################

"""Trait object which based on static dispatch determines how to recode from
the encoding of the source sequence to the encoding of the kmer"""
abstract type RecodingScheme end

"We can copy the encoding straight from the source to the kmer"
struct Copyable <: RecodingScheme end

"We can copy the encoding, then bitshift to create 4-bit encoding"
struct TwoToFour <: RecodingScheme end

"We skip all symbols in source that contain unmappable symbols"
struct FourToTwo <: RecodingScheme end

"We can use `BioSequences.ascii_encode`"
struct AsciiEncode <: RecodingScheme end

"The source is a bytevector, and we have no static knowledge of efficient
conversion to the right encoding"
struct GenericRecoding <: RecodingScheme end

function RecodingScheme(A::Alphabet, source_type::Type)::RecodingScheme
    return if source_type <: BioSequence
        if BioSequences.encoded_data_eltype(source_type) <: BitInteger
            As = Alphabet(source_type)
            if As == A
                Copyable()
            elseif As isa TwoBit && A isa TwoBit
                Copyable()
            elseif As isa FourBit && A isa FourBit
                Copyable()
            elseif As isa FourBit && A isa TwoBit
                FourToTwo()
            elseif As isa TwoBit && A isa FourBit
                TwoToFour()
            else
                GenericRecoding()
            end
        else
            GenericRecoding()
        end
    elseif source_type <: Bytes && BioSequences.codetype(A) isa BioSequences.AsciiAlphabet
        AsciiEncode()
    else
        GenericRecoding()
    end
end

################################################
# Unsafe extract
################################################

"Extract a full kmer at a given index of a sequence.
Note: These methods don't do any bounds checking"
function unsafe_extract end

@inline function unsafe_extract(
    ::TwoToFour,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    for i in 1:ksize(T)
        encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, i)))
        (_, data) = leftshift_carry(data, 4, encoding)
    end
    T(unsafe, data)
end

@inline function unsafe_extract(
    ::FourToTwo,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    for i in 1:ksize(T)
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))::UInt
        if count_ones(encoding) != 1
            throw(
                BioSequences.EncodeError(
                    Alphabet(T),
                    reinterpret(eltype(seq), encoding % UInt8),
                ),
            )
        end
        (_, data) = leftshift_carry(data, 2, trailing_zeros(encoding) % UInt)
    end
    T(unsafe, data)
end

@inline function unsafe_extract(
    ::Copyable,
    ::Type{T},
    seq::BioSequence,
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(seq))
    for i in from_index:(from_index + ksize(T) - 1)
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))::UInt
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end

@inline function unsafe_extract(
    ::AsciiEncode,
    ::Type{T},
    seq::AbstractVector{UInt8},
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(T))
    @inbounds for i in from_index:(from_index + ksize(T) - 1)
        byte = seq[i]
        encoding = BioSequences.ascii_encode(Alphabet(T), byte)
        if encoding > 0x7f
            throw(BioSequences.EncodeError(Alphabet(T), byte))
        end
        (_, data) = leftshift_carry(data, bps, encoding % UInt)
    end
    T(unsafe, data)
end

@inline function unsafe_extract(
    ::GenericRecoding,
    ::Type{T},
    seq,
    from_index,
) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(T))
    @inbounds for i in 1:ksize(T)
        symbol = convert(eltype(T), seq[i])
        encoding = UInt(BioSequences.encode(Alphabet(T), symbol))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    T(unsafe, data)
end

################################################
# Constructors with full parameterisation
################################################

function Kmer{A, K, N}(x) where {A, K, N}
    check_kmer(Kmer{A, K, N})
    build_kmer(RecodingScheme(A(), typeof(x)), Kmer{A, K, N}, x)
end

# BioSequences support indexing and fast length checks
@inline function build_kmer(R::RecodingScheme, ::Type{T}, s::BioSequence) where {T}
    length(s) == ksize(T) || error("Length of sequence must be K elements to build Kmer")
    unsafe_extract(R, T, s, 1)
end

# LongSequence with compatible alphabet: Extract whole coding elements
@inline function build_kmer(::Copyable, ::Type{T}, s::LongSequence) where {T}
    length(s) == ksize(T) || error("Length of sequence must be K elements to build Kmer")
    bps = BioSequences.BitsPerSymbol(Alphabet(T))
    data = ntuple(i -> BioSequences.reversebits(@inbounds(s.data[i]), bps), Val{nsize(T)}())
    (_, data) = rightshift_carry(data, bits_unused(T), zero(UInt))
    T(unsafe, data)
end

# TODO: LongSubSeq with compatible alphabet
# Note: LongSequence may be UInt64 whereas kmers use UInt32

# For UTF8-strings combined with an ASCII kmer alphabet, we convert to byte vector
@inline function build_kmer(
    R::AsciiEncode,
    ::Type{T},
    s::Union{String, SubString{String}},
) where {T}
    build_kmer(R, T, codeunits(s))
end

# For byte vectors, we can build a kmer iff the kmer alphabet is AsciiAlphabet
@inline function build_kmer(R::AsciiEncode, ::Type{T}, s::AbstractVector{UInt8}) where {T}
    length(s) == ksize(T) || error("Length of sequence must be K elements to build Kmer")
    unsafe_extract(R, T, s, 1)
end

# The generic fallback - dispatch on whether we can check length once
@inline function build_kmer(R::RecodingScheme, T::Type, s)
    build_kmer(Base.IteratorSize(typeof(s)), R, T, s)
end

@inline function build_kmer(::Base.SizeUnknown, ::RecodingScheme, T::Type, s)
    data = zero_tuple(T)
    A = Alphabet(T)
    bps = BioSequences.bits_per_symbol(A)
    i = 0
    for element in s
        i += 1
        i > ksize(T) && error("Length of sequence must be K elements to build Kmer")
        symbol = convert(eltype(A), element)
        carry = UInt(BioSequences.encode(A, symbol))
        (_, data) = leftshift_carry(data, bps, carry)
    end
    i == ksize(T) || error("Length of sequence must be K elements to build Kmer")
    T(unsafe, data)
end

@inline function build_kmer(
    ::Union{Base.HasLength, Base.HasShape},
    ::RecodingScheme,
    T::Type,
    s,
)
    length(s) == ksize(T) || error("Length of sequence must be K elements to build Kmer")
    data = zero_tuple(T)
    A = Alphabet(T)
    bps = BioSequences.bits_per_symbol(A)
    for element in s
        symbol = convert(eltype(A), element)
        carry = UInt(BioSequences.encode(A, symbol))
        (_, data) = leftshift_carry(data, bps, carry)
    end
    T(unsafe, data)
end

################################################
# Derived constructors
################################################

Kmer{A, K}(x) where {A, K} = derive_type(Kmer{A, K})(x)
Kmer{A1}(x::Kmer{A2, K, N}) where {A1, A2, K, N} = Kmer{A1, K, N}(x)

function kmer(::Val{K}, s::BioSequence{A}) where {A, K}
    K isa Int || error("K must be an Int")
    Kmer{A, K}(s)
end

################################################
# Construct other types from Kmers
################################################

# TODO: LongSequence
# TODO: String

################################################
# String literals
################################################

"""
    @mer_str -> Kmer

Construct a `Kmer` from the given string. The macro must be used with a flag
after the input string, e.g. `d` in `mer"TAG"d` or `a` in `mer"PCW"a`, signifying
the alphabet of the kmer.
The flags `d = DNAAlphabet{2}`, `r = RNAAlphabet{2}` and `a = AminoAcidAlphabet`
are recognized.

Because the macro is resolved and the kmer is created at parse time,
the macro is type stable, and may be used in high performance code.

# Examples
```jldoctest
julia> mer"UGCUA"r
RNA 5-mer:
UGCUA

julia> mer"YKVSTEDLLKKR"a
AminoAcid 12-mer:
YKVSTEDLLKKR

julia> mer"TATTAGCA"d
DNA 8-mer:
TATTAGCA
```
"""
macro mer_str(seq, flag)
    trimmed = BioSequences.remove_newlines(seq)
    ncu = ncodeunits(trimmed)
    # Unlike @dna_str, we default to 2-bit alphabets, because kmers
    # by convention are usually 2-bit only
    if flag == "dna" || flag == "d"
        Kmer{DNAAlphabet{2}, ncu}(trimmed)
    elseif flag == "rna" || flag == "r"
        Kmer{RNAAlphabet{2}, ncu}(trimmed)
    elseif flag == "aa" || flag == "a"
        Kmer{AminoAcidAlphabet, ncu}(trimmed)
    else
        error("Invalid type flag: '$(flag)'")
    end
end
