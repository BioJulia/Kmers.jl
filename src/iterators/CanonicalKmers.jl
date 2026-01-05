"""
    FwRvIterator{A <: NucleicAcidAlphabet, K, S}

Iterates 2-tuples of `(forward, reverse_complement)` of every kmer of type
`Kmer{A, K}` from the underlying sequence, in order.
`S` signifies the type of the underlying sequence.
This is typically more efficient than iterating over a `FwKmers` and
computing `reverse_complement` on every element.

See also: [`FwKmers`](@ref), [`CanonicalKmers`](@ref)

# Examples:
```jldoctest
julia> collect(FwRvIterator{DNAAlphabet{4}, 3}("AGCGT"))
3-element Vector{Tuple{Mer{3, DNAAlphabet{4}, 1}, Mer{3, DNAAlphabet{4}, 1}}}:
 (AGC, GCT)
 (GCG, CGC)
 (CGT, ACG)

julia> collect(FwRvIterator{DNAAlphabet{2}, 3}("AGNGT"))
ERROR: cannot encode 0x4e (Char 'N') in DNAAlphabet{2}
[...]
```
"""
struct FwRvIterator{A <: NucleicAcidAlphabet, K, S}
    seq::S

    function FwRvIterator{A, K, S}(seq::S) where {A, K, S}
        K isa Int || error("K must be an Int")
        K > 0 || error("K must be at least 1")
        return new{A, K, S}(seq)
    end
end

"`FwRvDNAIterator{K, S}`: Alias for `FwRvIterator{DNAAlphabet{2}, K, S}`"
const FwRvDNAIterator{K, S} = FwRvIterator{DNAAlphabet{2}, K, S}

"`FwRvRNAIterator{K, S}`: Alias for `FwRvIterator{RNAAlphabet{2}, K, S}`"
const FwRvRNAIterator{K, S} = FwRvIterator{RNAAlphabet{2}, K, S}

source_type(::Type{FwRvIterator{A, K, S}}) where {A, K, S} = S
kmertype(::Type{<:FwRvIterator{A, K}}) where {A, K} = derive_type(Kmer{A, K})
kmertype(it::FwRvIterator) = kmertype(typeof(it))
Base.eltype(T::Type{<:FwRvIterator{A, K}}) where {A, K} =
    Tuple{K, K} where {K <: kmertype(T)}

@inline function Base.length(it::FwRvIterator{A, K, S}) where {A, K, S}
    src = used_source(RecodingScheme(A(), S), it.seq)
    return max(0, length(src) - K + 1)
end

FwRvIterator{A, K}(s) where {A <: Alphabet, K} = FwRvIterator{A, K, typeof(s)}(s)

@inline function Base.iterate(it::FwRvIterator{A, K, S}, state...) where {A, K, S}
    return iterate_kmer(RecodingScheme(A(), S), it, state...)
end

# For the first kmer, we extract it, then reverse complement.
# When it's not done incrementally, it's faster to RC the whole
# kmer at once.
@inline function iterate_kmer(R::RecodingScheme, it::FwRvIterator{A, K}) where {A, K}
    length(it.seq) < K && return nothing
    fw = unsafe_extract(R, kmertype(it), it.seq, 1)
    rv = reverse_complement(fw)
    return ((fw, rv), (fw, rv, K + 1))
end

# Here, we need to convert to an abstractvector
@inline function iterate_kmer(
        R::AsciiEncode,
        it::FwRvIterator{A, K, S},
    ) where {A <: NucleicAcidAlphabet, K, S}
    src = used_source(RecodingScheme(A(), S), it.seq)
    Base.require_one_based_indexing(src)
    length(src) < K && return nothing
    fw = unsafe_extract(R, kmertype(it), src, 1)
    rv = reverse_complement(fw)
    return ((fw, rv), (fw, rv, K + 1))
end

@inline function iterate_kmer(
        ::GenericRecoding,
        it::FwRvIterator,
        state::Tuple{Kmer, Kmer, Int},
    )
    (fw, rv, i) = state
    i > length(it.seq) && return nothing
    symbol = convert(eltype(fw), @inbounds it.seq[i])
    fw = shift(fw, symbol)
    rv = shift_first(rv, complement(symbol))
    return ((fw, rv), (fw, rv, i + 1))
end

@inline function iterate_kmer(
        ::Copyable,
        it::FwRvIterator{<:TwoBit, K, <:BioSequence{<:TwoBit}},
        state::Tuple{Kmer, Kmer, Int},
    ) where {K}
    (fw, rv, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))
    fw = shift_encoding(fw, encoding)
    rv = shift_first_encoding(rv, encoding ⊻ 0x03)
    return ((fw, rv), (fw, rv, i + 1))
end

@inline function iterate_kmer(
        ::Copyable,
        it::FwRvIterator{<:FourBit, K, <:BioSequence{<:FourBit}},
        state::Tuple{Kmer, Kmer, Int},
    ) where {K}
    (fw, rv, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))
    fw = shift_encoding(fw, encoding)
    rc_encoding =
        reinterpret(UInt8, complement(reinterpret(eltype(rv), encoding % UInt8))) % UInt
    rv = shift_first_encoding(rv, rc_encoding)
    return ((fw, rv), (fw, rv, i + 1))
end

@inline function iterate_kmer(::TwoToFour, it::FwRvIterator, state::Tuple{Kmer, Kmer, Int})
    (fw, rv, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))
    fw = shift_encoding(fw, left_shift(UInt(1), encoding))
    rv = shift_first_encoding(rv, left_shift(UInt(1), encoding ⊻ 0x03))
    return ((fw, rv), (fw, rv, i + 1))
end

@inline function iterate_kmer(
        ::FourToTwo,
        it::FwRvIterator{A, K, <:BioSequence},
        state::Tuple{Kmer, Kmer, Int},
    ) where {A, K}
    (fw, rv, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))::UInt
    isone(count_ones(encoding)) || throw_uncertain(Alphabet(fw), eltype(it.seq), encoding)
    enc = trailing_zeros(encoding) % UInt
    fw = shift_encoding(fw, enc)
    rv = shift_first_encoding(rv, enc ⊻ 0x03)
    return ((fw, rv), (fw, rv, i + 1))
end

@inline function iterate_kmer(
        ::AsciiEncode,
        it::FwRvIterator{A},
        state::Tuple{Kmer, Kmer, Int},
    ) where {A}
    src = used_source(RecodingScheme(A(), source_type(typeof(it))), it.seq)
    Base.require_one_based_indexing(src)
    (fw, rv, i) = state
    i > length(src) && return nothing
    byte = @inbounds src[i]
    encoding = BioSequences.ascii_encode(A(), byte)
    if encoding > 0x7f
        throw(BioSequences.EncodeError(A(), repr(byte)))
    end
    # Hopefully this branch is eliminated at compile time...
    rc_encoding = if Alphabet(fw) isa FourBit
        reinterpret(UInt8, complement(reinterpret(DNA, encoding)))
    elseif Alphabet(fw) isa TwoBit
        encoding ⊻ 0x03
    else
        error(
            "Complementing encoding of a Nucleotide AsciiAlphabet which is neither 2 or 4 " *
                "bits have not been implemented yet.",
        )
    end
    fw = shift_encoding(fw, encoding % UInt)
    rv = shift_first_encoding(rv, rc_encoding % UInt)
    return ((fw, rv), (fw, rv, i + 1))
end

"""
    CanonicalKmers{A <: NucleicAcidAlphabet, K, S} <: AbstractKmerIterator{A, K}

Iterator of canonical nucleic acid kmers. The result of this iterator is equivalent
to calling `canonical` on each value of a `FwKmers` iterator, but may be more
efficient.

!!! note
    When counting small kmers, it may be more efficient to count `FwKmers`,
    then call `canonical` only once per unique kmer.

Can be constructed more conventiently with the constructors `CanonicalDNAMers{K}(s)`
`CanonicalRNAMers{K}(s)`

# Examples:
```jldoctest
julia> collect(CanonicalRNAMers{3}("AGCGA"))
3-element Vector{Kmer{RNAAlphabet{2}, 3, 1}}:
 AGC
 CGC
 CGA
```
"""
struct CanonicalKmers{A <: NucleicAcidAlphabet, K, S} <: AbstractKmerIterator{A, K}
    it::FwRvIterator{A, K, S}
end

source_type(::Type{CanonicalKmers{A, K, S}}) where {A, K, S} = S
@inline Base.length(it::CanonicalKmers) = length(it.it)

# Constructors
function CanonicalKmers{A, K}(s::S) where {S, A <: NucleicAcidAlphabet, K}
    return CanonicalKmers{A, K, S}(FwRvIterator{A, K}(s))
end
function CanonicalKmers{A, K, S}(s::S) where {S, A <: NucleicAcidAlphabet, K}
    return CanonicalKmers{A, K, S}(FwRvIterator{A, K}(s))
end

"`CanonicalDNAMers{K, S}`: Alias for `CanonicalKmers{DNAAlphabet{2}, K, S}`"
const CanonicalDNAMers{K, S} = CanonicalKmers{DNAAlphabet{2}, K, S}

"`CanonicalRNAMers{K, S}`: Alias for `CanonicalKmers{RNAAlphabet{2}, K, S}`"
const CanonicalRNAMers{K, S} = CanonicalKmers{RNAAlphabet{2}, K, S}

@inline function Base.iterate(it::CanonicalKmers{A, K, S}, state...) where {A, K, S}
    it = iterate(it.it, state...)
    isnothing(it) && return nothing
    ((fw, rv), state) = it
    return (fw < rv ? fw : rv, state)
end
