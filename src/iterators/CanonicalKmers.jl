"""
    CanonicalKmers{A <: NucleicAcidAlphabet, K, S}

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
    it::FwKmers{A, K, S}
end

source_type(::Type{CanonicalKmers{A, K, S}}) where {A, K, S} = S
@inline Base.length(it::CanonicalKmers) = length(it.it)

# Constructors
function CanonicalKmers{A, K}(s::S) where {S, A <: NucleicAcidAlphabet, K}
    CanonicalKmers{A, K, S}(FwKmers{A, K}(s))
end
function CanonicalKmers{A, K, S}(s::S) where {S, A <: NucleicAcidAlphabet, K}
    CanonicalKmers{A, K, S}(FwKmers{A, K}(s))
end

const CanonicalDNAMers{K, S} = CanonicalKmers{DNAAlphabet{2}, K, S}
const CanonicalRNAMers{K, S} = CanonicalKmers{RNAAlphabet{2}, K, S}

@inline function Base.iterate(it::CanonicalKmers{A, K, S}, state...) where {A, K, S}
    iterate_kmer(RecodingScheme(A(), S), it, state...)
end

# For the first kmer, we extract it, then reverse complement.
# When it's not done incrementally, it's faster to RC the whole
# kmer at once.
@inline function iterate_kmer(R::RecodingScheme, it::CanonicalKmers)
    length(it.it.seq) < ksize(eltype(it)) && return nothing
    fw = unsafe_extract(R, eltype(it), it.it.seq, 1)
    rv = reverse_complement(fw)
    (fw < rv ? fw : rv, (fw, rv, ksize(eltype(it)) + 1))
end

# Here, we need to convert to an abstractvector
@inline function iterate_kmer(
    R::AsciiEncode,
    it::CanonicalKmers{A, K, S},
) where {A <: NucleicAcidAlphabet, K, S <: Bytes}
    src = used_source(RecodingScheme(A(), S), it.it.seq)
    Base.require_one_based_indexing(src)
    length(src) < ksize(eltype(it)) && return nothing
    fw = unsafe_extract(R, eltype(it), src, 1)
    rv = reverse_complement(fw)
    (fw < rv ? fw : rv, (fw, rv, ksize(eltype(it)) + 1))
end

@inline function iterate_kmer(
    ::GenericRecoding,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    (fw, rv, i) = state
    i > length(it.it.seq) && return nothing
    symbol = convert(eltype(fw), @inbounds it.it.seq[i])
    fw = shift(fw, symbol)
    rv = shift_first(rv, complement(symbol))
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::Copyable,
    it::CanonicalKmers{<:TwoBit, K, <:BioSequence{<:TwoBit}},
    state::Tuple{Kmer, Kmer, Int},
) where {K}
    (fw, rv, i) = state
    i > length(it.it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.it.seq, i))
    fw = shift_encoding(fw, encoding)
    rv = shift_first_encoding(rv, encoding ⊻ 0x03)
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::Copyable,
    it::CanonicalKmers{<:FourBit, K, <:BioSequence{<:FourBit}},
    state::Tuple{Kmer, Kmer, Int},
) where {K}
    (fw, rv, i) = state
    i > length(it.it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.it.seq, i))
    fw = shift_encoding(fw, encoding)
    rc_encoding =
        reinterpret(UInt8, complement(reinterpret(eltype(rv), encoding % UInt8))) % UInt
    rv = shift_first_encoding(rv, rc_encoding)
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::TwoToFour,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    (fw, rv, i) = state
    i > length(it.it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.it.seq, i))
    fw = shift_encoding(fw, left_shift(UInt(1), encoding))
    rv = shift_first_encoding(rv, left_shift(UInt(1), encoding ⊻ 0x03))
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::FourToTwo,
    it::CanonicalKmers{A, K, <:BioSequence},
    state::Tuple{Kmer, Kmer, Int},
) where {A, K}
    (fw, rv, i) = state
    i > length(it.it.seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(it.it.seq, i))::UInt
    if count_ones(encoding) != 1
        throw(
            BioSequences.EncodeError(
                Alphabet(fw),
                reinterpret(eltype(it.it.seq), encoding % UInt8),
            ),
        )
    end
    enc = trailing_zeros(encoding) % UInt
    fw = shift_encoding(fw, enc)
    rv = shift_first_encoding(rv, enc ⊻ 0x03)
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::AsciiEncode,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    src = used_source(
        RecodingScheme(Alphabet(eltype(it)), source_type(typeof(it))),
        it.it.seq,
    )
    Base.require_one_based_indexing(src)
    (fw, rv, i) = state
    i > length(src) && return nothing
    byte = @inbounds src[i]
    encoding = BioSequences.ascii_encode(Alphabet(eltype(it)), byte)
    if encoding > 0x7f
        throw(BioSequences.EncodeError(Alphabet(eltype(it)), repr(byte)))
    end
    # Hopefully this branch is eliminated at compile time...
    rc_encoding = if Alphabet(fw) isa FourBit
        reinterpret(UInt8, complement(reinterpret(DNA, encoding)))
    elseif Alphabet(fw) isa TwoBit
        encoding ⊻ 0x03
    else
        error("Unreachable")
    end
    fw = shift_encoding(fw, encoding % UInt)
    rv = shift_first_encoding(rv, rc_encoding % UInt)
    (fw < rv ? fw : rv, (fw, rv, i + 1))
end
