@inline function BioSequences.extract_encoded_element(seq::Kmer{A}, i::Integer) where A
    T = typeof(seq)
    bps = BioSequences.bits_per_symbol(A()) % UInt
    index = div((i + n_unused(T) - 1) % UInt, per_word_capacity(T) % UInt) + 1
    offset = mod(((elements_in_head(T) - i) * bps) % UInt, 8 * sizeof(UInt))
    mask = UInt(1) << bps - 1
    right_shift(@inbounds(seq.data[index]), offset) & mask
end

# This is usually type unstable, but in user code, users may use constant-folded ranges,
# e.g. f(x) = x[2:4]. In this case, we need it to compile to very efficient code.
# Hence, it MUST use @inline
@inline function Base.getindex(kmer::Kmer{A}, range::AbstractRange{<:Integer}) where A
    @boundscheck checkbounds(kmer, range)
    K = length(range)
    N = n_coding_elements(Kmer{A, K})
    T = Kmer{A, K, N}
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for i in range
        (_, data) = leftshift_carry(data, nbits, BioSequences.extract_encoded_element(kmer, i))
    end
    T(unsafe, data)
end

# Same as above: This needs to be able to inline if the indices are known statically
@inline function Base.getindex(kmer::Kmer{A}, indices::AbstractVector{Bool}) where A
    @boundscheck checkbounds(eachindex(kmer), indices)
    K = sum(indices)
    N = n_coding_elements(Kmer{A, K})
    T = Kmer{A, K, N}
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for (i, bool) in enumerate(indices)
        bool || continue
        (_, data) = leftshift_carry(data, nbits, BioSequences.extract_encoded_element(kmer, i))
    end
    T(unsafe, data)
end

@inline function BioSequences.bitindex(kmer::Kmer, i::Unsigned)::Tuple{UInt, UInt}
    bps = BioSequences.bits_per_symbol(Alphabet(kmer)) % UInt
    bpe = (8 * sizeof(UInt)) % UInt
    (i, o) = divrem((UInt(i) - UInt(1) + n_unused(typeof(kmer))) * bps, bpe)
    o = bpe - o - bps
    i + 1, o
end

@inline function setindex(kmer::Kmer, i::Integer, s)
    @boundscheck checkbounds(kmer, i)
    bps = BioSequences.bits_per_symbol(Alphabet(kmer))
    symbol = convert(eltype(kmer), s)
    encoding = UInt(BioSequences.encode(Alphabet(kmer), symbol))
    (i, o) = BioSequences.bitindex(kmer, i % UInt)
    element = @inbounds kmer.data[i]
    mask = left_shift(UInt(1) << bps - 1, o)
    element &= ~mask
    element |= left_shift(encoding, o)
    typeof(kmer)(unsafe, @inbounds Base.setindex(kmer.data, element, i))
end