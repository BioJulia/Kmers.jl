################################################
# Trait dispatch
################################################

"""
    RecodingScheme

Trait object which determines the methods used to recode
data from one sequence into a `BioSequence`.
Any given construction of a kmer from a source sequence will
dispatch to one of the concrete subtypes of `RecodingScheme`
to determine how the data is copied most effectively.
"""
abstract type RecodingScheme end

"""
    Copyable <: RecodingScheme

Trait object that signifies that two sequences share identical encodings,
such that the encoded sequence can be copied directly.
This is the case for e.g. sequences of the same alphabet, or two
sequences of the alphabets `DNAAlphabet{2}` and `RNAAlphabet{2}`
"""
struct Copyable <: RecodingScheme end

"""
    TwoToFour <: RecodingScheme

Trait object that signifies a recoding from a 4-bit nucleotide sequence
to a 2-bit nucleotide sequence.
This can be more efficient than using `GenericRecoding`
"""
struct TwoToFour <: RecodingScheme end

"""
    FourToTwo <: RecodingScheme

Trait object that signifies a recoding from a 2-bit nucleotide sequence
to a 4-bit nucleotide sequence.
This can be more efficient than using `GenericRecoding`
"""
struct FourToTwo <: RecodingScheme end

"""
    AsciiEncode <: RecodingScheme

Trait object that signifies a recoding from an ASCII-encoded
`AbstractVector{UInt8}` to an `Alphabet` that implements the
`BioSequences.AsciiAlphabet` trait.
The validity of the bytes are checked during the recoding, including
checking that the sequence is actually ASCII encoded.
"""
struct AsciiEncode <: RecodingScheme end

"""
    GenericRecoding <: RecodingScheme

Trait object that signifies that none of the other recoding schemes
apply, that the fallback implementations must be used instead.
"""
struct GenericRecoding <: RecodingScheme end

"""
    is_ascii(::Type{T})::Bool

Trait function. Should return `true` for `AbstractVector{UInt8}`, or for
string types for which `codeunits(s)` returns an `AbstractVector{UInt8}`, where
every ASCII byte in the string is perserved in the vector.
This is true for all UTF8, latin1 and ASCII encoded string types.
"""
is_ascii(::Type) = false
is_ascii(::Type{<:Union{String, SubString{String}}}) = true
is_ascii(::Type{<:AbstractVector{UInt8}}) = true

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
    elseif is_ascii(source_type) && BioSequences.codetype(A) isa BioSequences.AsciiAlphabet
        AsciiEncode()
    else
        GenericRecoding()
    end
end

include("construction_utils.jl")

################################################
# Unsafe extract
################################################

@noinline function throw_uncertain(A::Alphabet, T::Type{<:BioSymbol}, enc::Unsigned)
    throw(BioSequences.EncodeError(A, reinterpret(T, enc % UInt8)))
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

################################################
# Construct other types from Kmers
################################################

function BioSequences.LongSequence{A}(kmer::Kmer{A}) where {A <: Alphabet}
    ratio = div(sizeof(UInt64), sizeof(UInt))
    data = zeros(UInt64, cld(n_coding_elements(typeof(kmer)), ratio))
    @inbounds for i in 1:n_coding_elements(typeof(kmer))
        target = cld(i, ratio)
        src = kmer.data[i] % UInt64
        if ratio == 2 && iseven(i)
            src <<= 32
        end
        data[target] |= src
    end
    bps = BioSequences.BitsPerSymbol(A())
    @inbounds for i in eachindex(data)
        data[i] = BioSequences.reversebits(data[i], bps)
    end
    bu = bits_unused(typeof(kmer))
    @inbounds if !iszero(bu)
        data[end] = data[end] >> bu
    end
    LongSequence{A}(data, length(kmer) % UInt)
end

# TODO: Do we want specialized constructors to contruct cross-alphabet longseqs
# from kmers?

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
