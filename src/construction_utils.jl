# This file contains construction utilities that are public, allowing users
# to create custom K-mer types such as syncmers, strobemers, minimizers etc.

# unsafe_extract
# shift_encoding
# shift_seq

"""
    unsafe_extract(::RecodingScheme, T::Type{<:Kmer}, seq, from::Int) -> T

Extract a `Kmer` of type `T` from `seq` beginning at index `from`.
This function is useful to create kmer or kmer-like types.

This function does not do any bounds checking, so the user must know
that `from:from+K-1` is inbounds in `seq`.
The validity of the data in the `seq` is validated by this function.

# Examples
```jldoctest
julia> seq = b"TAGCTAGA";

julia> Kmers.unsafe_extract(Kmers.AsciiEncode(), DNAKmer{4, 1}, seq, 2)
DNA 4-mer:
AGCT
```
"""
@inline function unsafe_extract(
        ::TwoToFour,
        ::Type{T},
        seq::BioSequence,
        from_index,
    ) where {T <: Kmer}
    data = zero_tuple(T)
    for i in from_index:(from_index + ksize(T) - 1)
        encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, i)))
        (_, data) = leftshift_carry(data, 4, encoding)
    end
    return T(unsafe, data)
end

@inline function unsafe_extract(
        ::FourToTwo,
        ::Type{T},
        seq::BioSequence,
        from_index,
    ) where {T <: Kmer}
    data = zero_tuple(T)
    for i in from_index:(from_index + ksize(T) - 1)
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))::UInt
        isone(count_ones(encoding)) || throw_uncertain(Alphabet(T), eltype(seq), encoding)
        (_, data) = leftshift_carry(data, 2, trailing_zeros(encoding) % UInt)
    end
    return T(unsafe, data)
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
    return T(unsafe, data)
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
    return T(unsafe, data)
end

@inline function unsafe_extract(
        ::GenericRecoding,
        ::Type{T},
        seq,
        from_index,
    ) where {T <: Kmer}
    data = zero_tuple(T)
    bps = BioSequences.bits_per_symbol(Alphabet(T))
    @inbounds for i in 0:(ksize(T) - 1)
        symbol = convert(eltype(T), seq[from_index + i])
        encoding = UInt(BioSequences.encode(Alphabet(T), symbol))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    return T(unsafe, data)
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
            throw_uncertain(Alphabet(kmer), eltype(seq), encoding)
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
