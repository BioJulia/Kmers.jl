struct CanonicalKmers{A <: Alphabet, K, S} <: AbstractKmerIterator{A, K}
    it::FwKmers{A, K, S}
end

source_type(::Type{CanonicalKmers{A, K, S}}) where {A, K, S} = S
load_source(x::CanonicalKmers) = x.it.seq
Base.length(it::CanonicalKmers) = length(it.it)

# Constructors
function CanonicalKmers{A, K}(s::S) where {S, A <: Alphabet, K}
    CanonicalKmers{S, A, K}(FwKmers{A, K}(s))
end

# Iteration
function Base.iterate(it::CanonicalKmers, state...)
    iterate_kmer(RecodingScheme(typeof(it)), it, state...)
end

# For these recoding schemes, no symbols in the source sequence are skipped.
# Hence, we can forward to just `extract`.
# Here, instead of reverse complementing each symbol, it's more efficient
# to do it in bulk by RC'ing the entire kmer
@inline function iterate_kmer(
    R::Union{GenericAlphabet, Copyable, TwoToFour, AsciiEncode, GenericBytes},
    it::CanonicalKmers,
)
    src = usable_source(it)
    length(src) < ksize(eltype(it)) && return nothing
    fw = extract(R, eltype(it), src, 1)
    rv = reverse_complement(fw)
    (min(fw, rv), (fw, rv, ksize(eltype(it)) + 1))
end

# Fallback: Just because it's Copyable doesn't mean we have neat bit-tricks
# to RC the encoding
@inline function iterate_kmer(
    ::Union{GenericAlphabet, Copyable},
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    src = usable_source(it)
    (fw, rv, i) = state
    i > length(src) && return nothing
    symbol = @inbounds src[i]
    encoding = UInt(BioSequences.encode(Alphabet(typeof(fw)), symbol))::UInt
    rc_encoding = UInt(BioSequences.encode(Alphabet(typeof(fw)), complement(symbol)))::UInt
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::Copyable,
    it::CanonicalKmers{<:TwoBit, K, <:TwoBit},
    state::Tuple{Kmer, Kmer, Int},
) where {K}
    src = usable_source(it)
    (fw, rv, i) = state
    i > length(src) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(src, i))::UInt
    rc_encoding = encoding ⊻ UInt(3)
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::Copyable,
    it::CanonicalKmers{<:FourBit, K, <:FourBit},
    state::Tuple{Kmer, Kmer, Int},
) where {K}
    src = usable_source(it)
    (fw, rv, i) = state
    i > length(src) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(src, i))::UInt
    rc_encoding = (@inbounds(FOURBIT_COMPLEMENT_LUT[encoding + UInt(1)])) % UInt
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::TwoToFour,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    src = usable_source(it)
    (fw, rv, i) = state
    i > length(src) && return nothing
    twobit_encoding = UInt(BioSequences.extract_encoded_element(src, i))::UInt
    fw_encoding = reinterpret(UInt8, decode(Alphabet(fw), twobit_encoding)) % UInt
    rc_encoding = reinterpret(UInt8, decode(Alphabet(fw), twobit_encoding ⊻ UInt(3))) % UInt
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

# 4 -> 2 (skipping): Ascii skipping LUT, as with FwKmers - DEFAULT STATE
@inline function iterate_kmer(
    ::Skipping,
    it::CanonicalKmers{A, K},
    state::Tuple{Kmer, Kmer, Int, Int}=(zero_kmer(Kmer{A, K}), zero_kmer(Kmer{A, K}), K, 1),
) where {A, K}
    src = usable_source(it)
    (fw, rv, remaining, i) = state
    while !iszero(remaining)
        i > length(src) && return nothing
        encoding = UInt(BioSequences.extract_encoded_element(src, i))::UInt
        i += 1
        if isone(count_ones(encoding))
            fw_encoding = trailing_zeros(encoding) % UInt
            fw = shift_encoding(fw, fw_encoding)
            rv = shift_first_encoding(rv, fw_encoding ⊻ UInt(3))
            remaining -= 1
        else
            remaining = K
            # No need to RC anything
            continue
        end
    end
    return (min(fw, rv), (fw, rv, 1, i))
end

@inline function iterate_kmer(
    ::AsciiSkipping,
    it::CanonicalKmers{A, K},
    state::Tuple{Kmer, Kmer, Int, Int}=(zero_kmer(Kmer{A, K}), zero_kmer(Kmer{A, K}), K, 1),
) where {A, K}
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (fw, rv, remaining, i) = state
    while !iszero(remaining)
        i > length(src) && return nothing
        byte = @inbounds src[i]
        i += 1
        encoding = @inbounds BYTE_LUT[byte + 0x01]
        encoding == 0xff && throw_bad_byte_error(byte)
        if encoding == 0xf0
            remaining = K
            continue
        else
            fw = shift_encoding(fw, encoding)
            rv = shift_first_encoding(rv, encoding ⊻ UInt(3))
            remaining -= 1
        end
    end
    return (min(fw, rv), (fw, rv, 1, i))
end

@inline function iterate_kmer(
    ::AsciiEncode,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (fw, rv, i) = state
    A = Alphabet(typeof(fw))
    i > length(src) && return nothing
    encoding = UInt(BioSequences.ascii_encode(A, @inbounds(src[i])))::UInt
    rc_encoding =
        UInt(BioSequences.encode(A, complement(BioSequences.decode(A, encoding))))::UInt
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::AsciiEncode,
    it::CanonicalKmers{<:FourBit},
    state::Tuple{Kmer, Kmer, Int},
)
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (fw, rv, i) = state
    A = Alphabet(typeof(fw))
    i > length(src) && return nothing
    encoding = UInt(BioSequences.ascii_encode(A, @inbounds(src[i])))::UInt
    rc_encoding = @inbounds(FOURBIT_COMPLEMENT_LUT[encoding + 0x01]) % UInt
    fw = shift_encoding(kmer, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rc), (fw, rv, i + 1))
end

@inline function iterate_kmer(
    ::GenericBytes,
    it::CanonicalKmers,
    state::Tuple{Kmer, Kmer, Int},
)
    src = usable_source(it)
    Base.require_one_based_indexing(src)
    (fw, rv, i) = state
    i > length(src) && return nothing
    char = reinterpret(Char, (src[i] % UInt32) << 24)
    fw_symbol = eltype(fw)(char)
    rc_symbol = complement(fw_symbol)
    fw = shift(fw, fw_symbol)
    rv = shift(rv, rc_symbol)
    (min(fw, rc), (fw, rv, i + 1))
end
