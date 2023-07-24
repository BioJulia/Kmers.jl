"""
    EveryKmer{S, A <: Alphabet, K}

Iterator of every forward kmer. `S` signifies the type of the underlying sequence,
and the eltype of the iterator is `Kmer{A, K, N}` with the appropriate `N`.

Can be constructed more conventiently with the constructors `EveryDNAMer{S, K}(s)`
and `EveryDNAMer{K}(s)`, and similar also for `EveryRNAMer` and `EveryAAMer`.

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
struct EveryKmer{S, A <: Alphabet, K} <: AbstractKmerIterator{A, K}
    seq::S
end

# Constructors
EveryKmer{A, K}(s) where {A <: Alphabet, K} = EveryKmer{typeof(s), A, K}
const EveryDNAMer{S, K} = EveryKmer{S, DNAAlphabet{2}, K}
const EveryRNAMer{S, K} = EveryKmer{S, RNAAlphabet{2}, K}
const EveryAAMer{S, K} = EveryKmer{S, AminoAcidAlphabet, K}

EveryDNAMer{K}(s) where K = EveryDNAMer{typeof(s), K}(s)
EveryRNAMer{K}(s) where K = EveryRNAMer{typeof(s), K}(s)
EveryAAMer{K}(s) where K = EveryAAMer{typeof(s), K}(s)

function EveryKmer{S, A, K}(s::S) where {S <: Union{String, SubString{String}}, A <: Alphabet, K}
    s2 = codeunits(s)
    EveryKmer{typeof(s2), A, K}(s2)
end

const SameEveryKmer{S, A, K} = EveryKmer{S, A} where {A, S <: BioSequence{A}}
const FourBit = Union{DNAAlphabet{4}, RNAAlphabet{4}}
const TwoBit = Union{DNAAlphabet{2}, RNAAlphabet{2}}

# Known length if every symbol of the sequence can be represented in the kmer
Base.IteratorSize(::Type{<:SameEveryKmer}) = Base.HasLength()
Base.IteratorSize(::Type{<:EveryKmer{<:BioSequence{<:TwoBit}, <:FourBit}}) = Base.HasLength()

function Base.length(it::SameEveryKmer{S, A, K}) where {S, A, K}
    length(it.seq) - K + 1
end

# These methods can carry the encoding directly over
function Base.iterate(it::EveryKmer{S, A, K}) where {A, K, S <: BioSequence{A}}
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

function Base.iterate(it::EveryKmer{S, A, K}, state::Tuple{Kmer, Integer}) where {A, K, S <: BioSequence{A}}
    seq = it.seq
    (kmer, i) = state
    i > length(seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(seq, i))
    new_kmer = q_push_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i+1))
end

# These methods can use special 2 -> 4 bit recoding
@inline recode(encoding::UInt) = left_shift(UInt(1), encoding)

function Base.iterate(it::EveryKmer{S, <:FourBit, K}) where {S <: BioSequence{<:TwoBit}, K}
    seq = it.seq
    length(seq) < K && return nothing
    data = zero_tuple(eltype(it))
    for i in 1:K
        encoding = recode(UInt(BioSequences.extract_encoded_element(seq, i)))
        (_, data) = leftshift_carry(data, 4, encoding)
    end
    kmer = eltype(it)(unsafe, data)
    (kmer, (kmer, K+1))
end

# TODO: Lots of code sharing in this file... can we refactor to be more clever?
function Base.iterate(it::EveryKmer{S, <:FourBit}, state::Tuple{Kmer, Integer}) where {S <: BioSequence{<:TwoBit}}
    seq = it.seq
    (kmer, i) = state
    i > length(seq) && return nothing
    encoding = recode(UInt(BioSequences.extract_encoded_element(seq, i)))
    new_kmer = q_push_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i+1))
end

# This is special because, by convention, we skip every ambiguous kmer
# instead of erroring.
function Base.iterate(
    it::EveryKmer{S, A, K}, state=(zero_kmer(Kmer{A, K}), K, 1)
) where {A <: TwoBit, S <: BioSequence{<:FourBit}, K}
    (kmer, remaining, i) = state
    seq = it.seq
    while !iszero(remaining)
        i > length(seq) && return nothing
        # TODO: Also, LUT here?
        encoding = UInt(BioSequences.extract_encoded_element(seq, i))
        i += 1
        # TODO: Is lookup table faster?
        remaining = ifelse(isone(count_ones(encoding)), remaining - 1, K)
        kmer = q_push_encoding(kmer, trailing_zeros(encoding) % UInt)
    end
    return (kmer, (kmer, 1, i))
end

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

# TODO: Change to lazy_str when new Julia LTS drops after 1.6
@noinline throw_bad_byte_error(b::UInt8) = error("Cannot interpret byte $(repr(b)) as nucleotide")

function Base.iterate(
    it::EveryKmer{S, A, K}, state=(zero_kmer(Kmer{A, K}), K, 1)
) where {A <: TwoBit, S <: AbstractVector{UInt8}, K}
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
        kmer = q_push_encoding(kmer, encoding % UInt)
    end
    return (kmer, (kmer, 1, i))
end