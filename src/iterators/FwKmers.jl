# TODO: Lots of code sharing in this file... can we refactor to be more clever?

"""
    FwKmers{A <: Alphabet, K, S}

Iterator of forward kmers. `S` signifies the type of the underlying sequence,
and the eltype of the iterator is `Kmer{A, K, N}` with the appropriate `N`.

Can be constructed more conventiently with the constructors `FwDNAMers{K}(s)`
and similar also for `FwRNAMers` and `FwAAMers`.

If `A <: Union{DNAAlphabet{2}, RNAAlphabet{2}}` and
`Alphabet(S) isa Union{DNAAlphabet{4}, RNAAlphabet{4}}`, the iterator skips all
kmers containing symbols not permitted in the 2-bit nucleotide alphabet.

# Examples:
```jldoctest
julia> v = collect(EveryDNAMer{3}("AGCGTATA"));

julia eltype(v), length(v)
(Kmer{DNAAlphabet{2}, 3, 1}, 6)

julia> length(collect(EveryRNAMer{3}(rna"UGDCUGAVC")))
2
```
"""
struct FwKmers{A <: Alphabet, K, S} <: AbstractKmerIterator{A, K}
    seq::S

    function FwKmers{A, K, S}(seq::S) where {A, K, S}
        K isa Int || error("K must be an Int")
        K > 0 || error("K must be at least 1")
        new{A, K, S}(seq)
    end
end

# Constructors
FwKmers{A, K}(s) where {A <: Alphabet, K} = FwKmers{A, K, typeof(s)}
const FwDNAMers{K, S} = FwKmers{DNAAlphabet{2}, K, S}
const FwRNAMers{K, S} = FwKmers{RNAAlphabet{2}, K, S}
const FwAAMers{K, S} = FwKmers{AminoAcidAlphabet, K, S}

FwDNAMers{K}(s) where K = FwDNAMers{K, typeof(s), }(s)
FwRNAMers{K}(s) where K = FwRNAMers{K, typeof(s)}(s)
FwAAMers{K}(s) where K = FwAAMers{K, typeof(s)}(s)

function FwKmers{A, K}(s::S) where {S <: Union{String, SubString{String}}, A <: Alphabet, K}
    s2 = codeunits(s)
    FwKmers{A, K, typeof(s2)}(s2)
end

# Known length if every symbol of the sequence can be represented in the kmer
Base.IteratorSize(::Type{<:FwKmers{A, K, <:BioSequence{A}}}) where {A <: Alphabet, K} = Base.HasLength()
Base.IteratorSize(::Type{<:FwKmers{<:FourBit, K, <:BioSequence{<:TwoBit}}}) where K = Base.HasLength()

function Base.length(it::FwKmers{A, K, <:BioSequence{A}}) where {A <: Alphabet, K}
    length(it.seq) - K + 1
end

# Generic fallback
function Base.iterate(it::FwKmers{A, K, S}) where {A <: Alphabet, K, S}
    seq = it.seq
    length(seq) < K && return nothing
    data = zero_tuple(eltype(it))
    bps = BioSequences.bits_per_symbol(A())
    @inbounds for i in 1:K
        symbol = seq[i]
        encoding = UInt(BioSequences.encode(A(), convert(eltype(A), symbol)))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    kmer = eltype(it)(unsafe, data)
    (kmer, (kmer, K+1))
end

function Base.iterate(it::FwKmers, state::Tuple{Kmer, Integer})
    seq = it.seq
    (kmer, i) = state
    i > length(seq) && return nothing
    symbol = @inbounds seq[i]
    new_kmer = shift(kmer, convert(eltype(A), symbol))
    (new_kmer, (new_kmer, i+1))
end

# These methods can carry the encoding directly over. We call into the internal method
# `iterate_copy`, because specifying the precise type constrains (either the same alphabet
# in the sequence and the iterator, OR both have either TwoBit or FourBit)
# is quite hard.
function Base.iterate(it::FwKmers{A, K, <:BioSequence{A}, }) where {A <: Alphabet, K}
    iterate_copy(it)
end

function Base.iterate(it::FwKmers{<:TwoBit, K, <:BioSequence{<:TwoBit}}) where K
    iterate_copy(it)
end

function Base.iterate(it::FwKmers{<:FourBit, K, <:BioSequence{<:FourBit}}) where K
    iterate_copy(it)
end

function Base.iterate(it::FwKmers{A, K, <:BioSequence{A}}, state::Tuple{Kmer, Integer}) where {A <: Alphabet, K}
    iterate_copy(it, state)
end

function Base.iterate(it::FwKmers{<:TwoBit, K, <:BioSequence{<:TwoBit}}, state::Tuple{Kmer, Integer}) where K
    iterate_copy(it, state)
end

function Base.iterate(it::FwKmers{<:FourBit, K, <:BioSequence{<:FourBit}}, state::Tuple{Kmer, Integer}) where K
    iterate_copy(it, state)
end

@inline function iterate_copy(it::FwKmers{A, K, S}) where {A, K, S}
    seq = it.seq
    length(seq) < K && return nothing
    data = zero_tuple(eltype(it))
    bps = BioSequences.bits_per_symbol(A())
    for i in 1:K
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    kmer = eltype(it)(unsafe, data)
    (kmer, (kmer, K+1))
end

@inline function iterate_copy(it::FwKmers, state::Tuple{Kmer, Integer})
    seq = it.seq
    (kmer, i) = state
    i > length(seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(seq, i))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i+1))
end

# These methods can use special 2 -> 4 bit recoding
function Base.iterate(it::FwKmers{<:FourBit, K, S}) where {S <: BioSequence{<:TwoBit}, K}
    seq = it.seq
    length(seq) < K && return nothing
    data = zero_tuple(eltype(it))
    for i in 1:K
        encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, i)))
        (_, data) = leftshift_carry(data, 4, encoding)
    end
    kmer = eltype(it)(unsafe, data)
    (kmer, (kmer, K+1))
end

function Base.iterate(it::FwKmers{<:FourBit, K, S}, state::Tuple{Kmer, Integer}) where {K, S <: BioSequence{<:TwoBit}}
    seq = it.seq
    (kmer, i) = state
    i > length(seq) && return nothing
    encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(seq, i)))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i+1))
end

# This is special because, by convention, we skip every ambiguous kmer
# instead of erroring.
function Base.iterate(
    it::FwKmers{A, K, S}, state=(zero_kmer(Kmer{A, K}), K, 1)
) where {A <: TwoBit, K, S <: BioSequence{<:FourBit}}
    (kmer, remaining, i) = state
    seq = it.seq
    while !iszero(remaining)
        i > length(seq) && return nothing
        # TODO: Also, LUT here?
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))
        i += 1
        # TODO: Is lookup table faster?
        remaining = ifelse(isone(count_ones(encoding)), remaining - 1, K)
        kmer = shift_encoding(kmer, trailing_zeros(encoding) % UInt)
    end
    return (kmer, (kmer, 1, i))
end

function Base.iterate(
    it::FwKmers{A, K}, state=(zero_kmer(Kmer{A, K}), K, 1)
) where {A <: TwoBit, K}
    (kmer, remaining, i) = state
    seq = it.seq
    Base.require_one_based_indexing(seq)
    while !iszero(remaining)
        i > length(seq) && return nothing
        byte = @inbounds seq[i]
        i += 1
        encoding = @inbounds BYTE_LUT[byte + 0x01]
        encoding == 0xff && throw_bad_byte_error(byte)
        remaining = ifelse(encoding == 0xf0, K, remaining - 1)
        kmer = shift_encoding(kmer, encoding % UInt)
    end
    return (kmer, (kmer, 1, i))
end
