"""
    FwKmers{A <: Alphabet, K, S}

Iterator of forward kmers. `S` signifies the type of the underlying sequence,
and the eltype of the iterator is `Kmer{A, K, N}` with the appropriate `N`.

Can be constructed more conventiently with the constructors `FwDNAMers{K}(s)`
and similar also for `FwRNAMers` and `FwAAMers`.

# Examples:
```jldoctest
julia> v = collect(FwDNAMers{3}("AGCGTATA"));

julia eltype(v), length(v)
(Kmer{DNAAlphabet{2}, 3, 1}, 6)

julia> length(collect(FwRNAMers{3}(rna"UGDCUGAVC")))
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

source_type(::Type{FwKmers{A, K, S}}) where {A, K, S} = S
load_source(x::FwKmers) = x.seq

function Base.length(it::FwKmers{A, K, S}) where {A, K, S}
    src = used_source(RecodingScheme(A(), S), it.seq)
    length(src) - ksize(eltype(it)) + 1
end

# Constructors
FwKmers{A, K}(s) where {A <: Alphabet, K} = FwKmers{A, K, typeof(s)}

const FwDNAMers{K, S} = FwKmers{DNAAlphabet{2}, K, S}
const FwRNAMers{K, S} = FwKmers{RNAAlphabet{2}, K, S}
const FwAAMers{K, S} = FwKmers{AminoAcidAlphabet, K, S}

FwDNAMers{K}(s) where {K} = FwDNAMers{K, typeof(s)}(s)
FwRNAMers{K}(s) where {K} = FwRNAMers{K, typeof(s)}(s)
FwAAMers{K}(s) where {K} = FwAAMers{K, typeof(s)}(s)

function Base.iterate(it::FwKmers{A, K, S}, state...) where {A, K, S}
    iterate_kmer(RecodingScheme(A(), S), it, state...)
end

# For the first kmer, we just forward to `unsafe_extract`
@inline function iterate_kmer(R::RecodingScheme, it::FwKmers)
    length(it.seq) < ksize(eltype(it)) && return nothing
    kmer = unsafe_extract(R, eltype(it), it.seq, 1)
    (kmer, (kmer, ksize(eltype(it)) + 1))
end

# Here, we need to convert to an abstractvector
@inline function iterate_kmer(
    R::AsciiEncode,
    it::FwKmers{A, K, S},
) where {A <: Alphabet, K, S <: Bytes}
    src = used_source(RecodingScheme(A(), S), it.seq)
    Base.require_one_based_indexing(src)
    length(src) < ksize(eltype(it)) && return nothing
    kmer = unsafe_extract(R, eltype(it), src, 1)
    (kmer, (kmer, ksize(eltype(it)) + 1))
end

@inline function iterate_kmer(::GenericRecoding, it::FwKmers, state::Tuple{Kmer, Int})
    (kmer, i) = state
    i > length(it.seq) && return nothing
    symbol = @inbounds it.seq[i]
    new_kmer = shift(kmer, convert(eltype(kmer), symbol))
    (new_kmer, (new_kmer, nextind(it.seq, i)))
end

@inline function iterate_kmer(::Copyable, it::FwKmers, state::Tuple{Kmer, Int})
    (kmer, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, nextind(it.seq, i)))
end

@inline function iterate_kmer(::TwoToFour, it::FwKmers, state::Tuple{Kmer, Int})
    (kmer, i) = state
    i > length(it.seq) && return nothing
    encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(it.seq, i)))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, nextind(it.seq, i)))
end

@inline function iterate_kmer(
    ::FourToTwo,
    it::FwKmers{A, K, <:BioSequence},
    state::Tuple{Kmer, Int},
) where {A, K}
    (kmer, i) = state
    i > length(it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))::UInt
    if count_ones(encoding) != 1
        throw(
            BioSequences.EncodeError(
                Alphabet(kmer),
                reinterpret(eltype(it.seq), encoding % UInt8),
            ),
        )
    end
    kmer = shift_encoding(kmer, trailing_zeros(encoding) % UInt)
    return (kmer, (kmer, nextind(it.seq, i)))
end

@inline function iterate_kmer(::AsciiEncode, it::FwKmers, state::Tuple{Kmer, Int})
    src = used_source(RecodingScheme(Alphabet(eltype(it)), source_type(typeof(it))), it.seq)
    Base.require_one_based_indexing(src)
    (kmer, i) = state
    i > length(src) && return nothing
    byte = @inbounds src[i]
    encoding = BioSequences.ascii_encode(Alphabet(eltype(it)), byte)
    if encoding > 0x7f
        throw(BioSequences.EncodeError(Alphabet(eltype(it)), repr(byte)))
    end
    kmer = shift_encoding(kmer, encoding % UInt)
    return (kmer, (kmer, nextind(src, i)))
end
