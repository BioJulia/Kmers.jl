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

function Base.length(it::FwKmers)
    Base.IteratorSize(typeof(it)) == Base.HasLength() || throw(MethodError(length, (it,)))
    length(usable_source(it)) - ksize(eltype(it)) + 1
end

# Constructors
const FwDNAMers{K, S} = FwKmers{DNAAlphabet{2}, K, S}
const FwRNAMers{K, S} = FwKmers{RNAAlphabet{2}, K, S}
const FwAAMers{K, S} = FwKmers{AminoAcidAlphabet, K, S}

FwKmers{A, K}(s) where {A <: Alphabet, K} = FwKmers{A, K, typeof(s)}

FwDNAMers{K}(s) where {K} = FwDNAMers{K, typeof(s)}(s)
FwRNAMers{K}(s) where {K} = FwRNAMers{K, typeof(s)}(s)
FwAAMers{K}(s) where {K} = FwAAMers{K, typeof(s)}(s)

function Base.iterate(it::FwKmers, state...)
    iterate_kmer(RecodingScheme(typeof(it)), it, state...)
end

# For these recoding schemes, no symbols in the source sequence are skipped.
# Hence, we can forward to just `extract`
@inline function iterate_kmer(
    R::Union{GenericAlphabet, Copyable, TwoToFour, AsciiEncode, GenericBytes},
    it::FwKmers,
)
    src = usable_source(it)
    length(src) < ksize(eltype(it)) && return nothing
    kmer = extract(R, eltype(it), src, 1)
    (kmer, (kmer, ksize(eltype(it)) + 1))
end

@inline function iterate_kmer(::GenericAlphabet, it::FwKmers, state::Tuple{Kmer, Int})
    src = usable_source(it)
    (kmer, i) = state
    i > length(src) && return nothing
    symbol = @inbounds src[i]
    new_kmer = shift(kmer, convert(eltype(kmer), symbol))
    (new_kmer, (new_kmer, i + 1))
end

@inline function iterate_kmer(::Copyable, it::FwKmers, state::Tuple{Kmer, Int})
    src = usable_source(it)
    (kmer, i) = state
    i > length(src) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(src, i))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i + 1))
end

@inline function iterate_kmer(::TwoToFour, it::FwKmers, state::Tuple{Kmer, Int})
    src = usable_source(it)
    (kmer, i) = state
    i > length(src) && return nothing
    encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(src, i)))
    new_kmer = shift_encoding(kmer, encoding)
    (new_kmer, (new_kmer, i + 1))
end

@inline function iterate_kmer(
    ::Skipping,
    it::FwKmers{A, K},
    state::Tuple{Kmer, Int, Int}=(zero_kmer(Kmer{A, K}), K, 1),
) where {A, K}
    (kmer, remaining, i) = state
    src = usable_source(it)
    while !iszero(remaining)
        i > length(src) && return nothing
        encoding = UInt(BioSequences.extract_encoded_element(src, i))
        i += 1
        remaining = ifelse(isone(count_ones(encoding)), remaining - 1, K)
        kmer = shift_encoding(kmer, trailing_zeros(encoding) % UInt)
    end
    return (kmer, (kmer, 1, i))
end

@inline function iterate_kmer(::GenericBytes, it::FwKmers, state::Tuple{Kmer, Int})
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (kmer, i) = state
    i > length(src) && return nothing
    char = reinterpret(Char, (src[i] % UInt32) << 24)
    symbol = eltype(eltype(it))(char)
    kmer = shift(kmer, symbol)
    return (kmer, (kmer, i + 1))
end

@inline function iterate_kmer(::AsciiEncode, it::FwKmers, state::Tuple{Kmer, Int})
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (kmer, i) = state
    i > length(src) && return nothing
    encoding = BioSequences.ascii_encode(Alphabet(eltype(it)), @inbounds(src[i]))
    kmer = shift_encoding(kmer, encoding)
    return (kmer, (kmer, i + 1))
end

@inline function iterate_kmer(
    ::AsciiSkipping,
    it::FwKmers{A, K},
    state=(zero_kmer(Kmer{A, K}), K, 1),
) where {A, K}
    (kmer, remaining, i) = state
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    while !iszero(remaining)
        i > length(src) && return nothing
        byte = @inbounds src[i]
        i += 1
        encoding = @inbounds ASCII_SKIPPING_LUT[byte + 0x01]
        encoding == 0xff && throw_bad_byte_error(byte)
        remaining = ifelse(encoding == 0xf0, K, remaining - 1)
        kmer = shift_encoding(kmer, encoding % UInt)
    end
    return (kmer, (kmer, 1, i))
end
