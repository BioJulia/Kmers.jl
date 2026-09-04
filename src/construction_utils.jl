# This file contains construction utilities that are public, allowing users
# to create custom K-mer types such as syncmers, strobemers, minimizers etc.

# unsafe_extract
# shift_encoding
# shift_seq

################################################
# Per-symbol recoding
################################################

# Keep recoding of individual symbols in these small helpers so both Kmer and
# Oligo construction use exactly the same validation and encoding rules.
@inline function recode_element(
        ::TwoToFour,
        ::Alphabet,
        source::BioSequence,
        index::Integer,
        ::Type{U},
    ) where {U <: Unsigned}
    encoding = BioSequences.extract_encoded_element(source, index) % U
    return left_shift(one(U), encoding)
end

@inline function recode_element(
        ::FourToTwo,
        alphabet::Alphabet,
        source::BioSequence,
        index::Integer,
        ::Type{U},
    ) where {U <: Unsigned}
    encoding = BioSequences.extract_encoded_element(source, index) % U
    isone(count_ones(encoding)) || throw_uncertain(alphabet, eltype(source), encoding)
    return trailing_zeros(encoding) % U
end

@inline function recode_element(
        ::Copyable,
        ::Alphabet,
        source::BioSequence,
        index::Integer,
        ::Type{U},
    ) where {U <: Unsigned}
    return BioSequences.extract_encoded_element(source, index) % U
end

@inline function recode_element(
        ::AsciiEncode,
        alphabet::Alphabet,
        source::AbstractVector{UInt8},
        index::Integer,
        ::Type{U},
    ) where {U <: Unsigned}
    byte = @inbounds source[index]
    encoding = BioSequences.ascii_encode(alphabet, byte)
    if encoding > 0x7f
        throw(BioSequences.EncodeError(alphabet, byte))
    end
    return encoding % U
end

@inline function recode_element(
        ::GenericRecoding,
        alphabet::Alphabet,
        source,
        index::Integer,
        ::Type{U},
    ) where {U <: Unsigned}
    symbol = convert(eltype(alphabet), @inbounds(source[index]))
    return BioSequences.encode(alphabet, symbol) % U
end

@inline function extract_elements(
        recoding::RecodingScheme,
        ::Type{T},
        seq,
        from_index,
    ) where {T <: Kmer}
    data = zero_tuple(T)
    alphabet = Alphabet(T)
    bps = BioSequences.bits_per_symbol(alphabet)
    for offset in 0:(ksize(T) - 1)
        encoding = recode_element(recoding, alphabet, seq, from_index + offset, UInt)
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    return T(unsafe, data)
end

"""
    unsafe_extract(::RecodingScheme, T::Type{<:Kmer}, seq, from::Int) -> T
    unsafe_extract(::RecodingScheme, T::Type{<:Oligo}, seq, from::Int, len::Int) -> T

Extract a `Kmer` of type `T` from `seq` beginning at index `from`, or an
`Oligo` containing `len` symbols beginning at that index. This function is
useful to create kmer or kmer-like types.

This function does not do any bounds checking, so the user must know
that the extracted range is inbounds in `seq`. For `Oligo`, the user must
also know that `len` does not exceed [`capacity(T)`](@ref).
The validity of the data in the `seq` is validated by this function.

# Examples
```jldoctest
julia> seq = b"TAGCTAGA";

julia> Kmers.unsafe_extract(Kmers.AsciiEncode(), DNAKmer{4, 1}, seq, 2)
DNA 4-mer:
AGCT

julia> Kmers.unsafe_extract(Kmers.AsciiEncode(), DNAOligo{UInt16}, seq, 2, 4)
4nt DNAOligo{UInt16}:
AGCT
```
"""
@inline function unsafe_extract(
        recoding::RecodingScheme,
        ::Type{T},
        seq,
        from_index,
    ) where {T <: Kmer}
    return extract_elements(recoding, T, seq, from_index)
end

################################################
# Packed extraction from LongSequence and LongSubSeq
################################################

# Load `len` adjacent encoded symbols into the low bits of `U`, with the first
# symbol in the least significant field. Long sequences are backed by UInt64s,
# so wide backing integers are assembled from as many source words as needed.
# An empty extraction never reads `source.data`.
@inline function load_packed_bits(
        ::Type{UInt64},
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
        ::BioSequences.BitsPerSymbol{B},
    ) where {B}
    iszero(len) && return zero(UInt64)

    bit_index = BioSequences.bitindex(source, from)
    word_index = BioSequences.index(bit_index)
    bit_offset = BioSequences.offset(bit_index)
    nbits = len * B
    bits = right_shift(@inbounds(source.data[word_index]), bit_offset)
    if bit_offset + nbits > 64
        bits |= left_shift(@inbounds(source.data[word_index + 1]), 64 - bit_offset)
    end
    return bits & BioSequences.bitmask(UInt64, nbits)
end

@inline function load_packed_bits(
        ::Type{U},
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
        ::BioSequences.BitsPerSymbol{B},
    ) where {U <: Unsigned, B}
    iszero(len) && return zero(U)

    bit_index = BioSequences.bitindex(source, from)
    word_index = BioSequences.index(bit_index)
    source_offset = Int(BioSequences.offset(bit_index))
    remaining = len * B
    destination_offset = 0
    bits = zero(U)
    while remaining > 0
        nbits = min(remaining, 64 - source_offset)
        word = right_shift(@inbounds(source.data[word_index]), source_offset)
        bits |= left_shift((word & BioSequences.bitmask(UInt64, nbits)) % U, destination_offset)
        remaining -= nbits
        destination_offset += nbits
        word_index += 1
        source_offset = 0
    end
    return bits
end

@inline function packed_word(
        ::Copyable,
        ::Alphabet,
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
    )
    return load_packed_bits(UInt64, source, from, len, BioSequences.BitsPerSymbol(source))
end

@inline function packed_word(
        ::TwoToFour,
        ::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{2}},
        from::Int,
        len::Int,
    )
    bits = load_packed_bits(UInt64, source, from, len, BioSequences.BitsPerSymbol{2}())
    encoding = BioSequences.two_to_four_bits(bits % UInt32)
    return encoding & BioSequences.bitmask(UInt64, 4 * len)
end

@inline function check_four_to_two_word(
        alphabet::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        bits::UInt64,
        len::Int,
    )
    mask = BioSequences.bitmask(UInt64, 4 * len) & 0x1111111111111111
    uncertain = BioSequences.uncertain_kernel(Alphabet(source), bits) & mask
    if !iszero(uncertain)
        offset = trailing_zeros(uncertain)
        encoding = right_shift(bits, offset) & 0x0f
        throw_uncertain(alphabet, eltype(source), encoding)
    end
    return nothing
end

@inline function four_to_two_half(
        alphabet::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        from::Int,
        len::Int,
    )
    bits = load_packed_bits(UInt64, source, from, len, BioSequences.BitsPerSymbol{4}())
    check_four_to_two_word(alphabet, source, bits, len)
    return BioSequences.four_to_two_bits(bits) % UInt64
end

@inline function packed_word(
        ::FourToTwo,
        alphabet::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        from::Int,
        len::Int,
    )
    first_len = min(len, 16)
    first_word = four_to_two_half(alphabet, source, from, first_len)
    if len == first_len
        return first_word & BioSequences.bitmask(UInt64, 2 * len)
    end
    second_word = four_to_two_half(alphabet, source, from + first_len, len - first_len)
    return first_word | left_shift(second_word, 32)
end

@inline function extract_packed(
        recoding::Union{Copyable, TwoToFour, FourToTwo},
        ::Type{T},
        source::BioSequences.SeqOrView,
        from::Int,
        bps::BioSequences.BitsPerSymbol{B},
    ) where {T <: Kmer, B}
    iszero(B) && return T(unsafe, zero_tuple(T))
    nsize(T) == 0 && return T(unsafe, zero_tuple(T))

    alphabet = Alphabet(T)
    per_word = div(64, B)
    words = ntuple(Val(nsize(T))) do i
        offset = (i - 1) * per_word
        len = min(per_word, ksize(T) - offset)
        bits = packed_word(recoding, alphabet, source, from + offset, len)
        UInt(B == 1 ? Base.bitreverse(bits) : BioSequences.reversebits(bits, bps))
    end
    (_, data) = rightshift_carry(words, bits_unused(T), zero(UInt))
    return T(unsafe, data)
end

@inline function unsafe_extract(
        recoding::Copyable,
        ::Type{T},
        source::BioSequences.SeqOrView,
        from::Int,
    ) where {T <: Kmer}
    bps = BioSequences.BitsPerSymbol(Alphabet(T))
    if bps == BioSequences.BitsPerSymbol(source) && !iszero(BioSequences.bits_per_symbol(bps))
        return extract_packed(recoding, T, source, from, bps)
    end
    return extract_elements(recoding, T, source, from)
end

@inline function unsafe_extract(
        recoding::TwoToFour,
        ::Type{T},
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{2}},
        from::Int,
    ) where {T <: Kmer{<:NucleicAcidAlphabet{4}}}
    return extract_packed(recoding, T, source, from, BioSequences.BitsPerSymbol(Alphabet(T)))
end

@inline function unsafe_extract(
        recoding::FourToTwo,
        ::Type{T},
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        from::Int,
    ) where {T <: Kmer{<:NucleicAcidAlphabet{2}}}
    return extract_packed(recoding, T, source, from, BioSequences.BitsPerSymbol(Alphabet(T)))
end

##########################
# Shift encoding
##########################

"""
    shift_encoding(kmer::T, encoding::UInt) where {T <: Kmer} -> T

Add `encoding`, a valid encoding in the alphabet of the `kmer`,
to the end of `kmer` and discarding the first symbol in `kmer`.

It is the user's responsibility to ensure that `encoding` is valid.

# Examples
```jldoctest
julia> enc = UInt(0x0a); # encoding of DNA_Y in 4-bit alphabets

julia> kmer = Kmer{DNAAlphabet{4}, 4}("TAGA");

julia> Kmers.shift_encoding(kmer, enc)
DNA 4-mer:
AGAY
```
"""
@inline function shift_encoding(kmer::Kmer, encoding::UInt)
    isempty(kmer) && return kmer
    bps = BioSequences.bits_per_symbol(kmer)
    (_, new_data) = leftshift_carry(kmer.data, bps, encoding)
    return typeof(kmer)(unsafe, (first(new_data) & get_mask(typeof(kmer)), Base.tail(new_data)...))
end

###########################

"""
    unsafe_shift_from(::RecodingScheme, kmer::T, seq, from::Int, ::Val{S}) -> T

Extract `S::Int` symbols from sequence `seq` at positions `from:from+S-1`,
and shift them into `kmer`.

This function does not do any bounds checking, so it is the user's
responsibility to ensure that `from` is inbounds, and the recoding scheme
valid.
It is assumed that `S < K`, where `K == length(kmer)`. If `S ≥ K`, use
[`unsafe_extract`](@ref) instead.

# Examples
```jldoctest
julia> seq = dna"TAGCGGA";

julia> kmer = mer"GGTG"d;

julia> Kmers.unsafe_shift_from(Kmers.FourToTwo(), kmer, seq, 3, Val(2))
DNA 4-mer:
TGGC
```
"""
@inline function unsafe_shift_from(
        ::GenericRecoding,
        kmer::Kmer,
        seq,
        from::Int,
        ::Val{S},
    ) where {S}
    for i in 0:(S - 1)
        symbol = @inbounds seq[from + i]
        kmer = shift(kmer, convert(eltype(kmer), symbol))
    end
    return kmer
end

@inline function unsafe_shift_from(
        ::Copyable,
        kmer::Kmer,
        seq::BioSequence,
        from::Int,
        ::Val{S},
    ) where {S}
    for i in 0:(S - 1)
        encoding = UInt(BioSequences.extract_encoded_element(seq, from + i))
        kmer = shift_encoding(kmer, encoding)
    end
    return kmer
end

@inline function unsafe_shift_from(
        ::TwoToFour,
        kmer::Kmer{<:NucleicAcidAlphabet{4}},
        seq::BioSequence{<:NucleicAcidAlphabet{2}},
        from::Int,
        ::Val{S},
    ) where {S}
    for i in 0:(S - 1)
        encoding =
            left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, from + i)))
        kmer = shift_encoding(kmer, encoding)
    end
    return kmer
end

@inline function unsafe_shift_from(
        ::FourToTwo,
        kmer::Kmer{<:NucleicAcidAlphabet{2}},
        seq::BioSequence{<:NucleicAcidAlphabet{4}},
        from::Int,
        ::Val{S},
    ) where {S}
    for i in 0:(S - 1)
        encoding = UInt(BioSequences.extract_encoded_element(seq, from + i))::UInt
        isone(count_ones(encoding)) ||
            throw_uncertain(Alphabet(kmer)::Alphabet, eltype(seq)::Type{<:BioSymbol}, encoding)
        kmer = shift_encoding(kmer, trailing_zeros(encoding) % UInt)
    end
    return kmer
end

@inline function unsafe_shift_from(
        ::AsciiEncode,
        kmer::Kmer,
        seq::AbstractVector{UInt8},
        from::Int,
        ::Val{S},
    ) where {S}
    for i in 0:(S - 1)
        byte = @inbounds seq[from + i]
        encoding = BioSequences.ascii_encode(Alphabet(typeof(kmer)), byte)
        if encoding > 0x7f
            throw(BioSequences.EncodeError(Alphabet(typeof(kmer)), byte))
        end
        kmer = shift_encoding(kmer, encoding % UInt)
    end
    return kmer
end
