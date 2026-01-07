@inline function BioSequences.extract_encoded_element(seq::Kmer, i::Integer)
    T = typeof(seq)
    bps = BioSequences.bits_per_symbol(Alphabet(seq)) % UInt
    index = div((i + n_unused(T) - 1) % UInt, per_word_capacity(T) % UInt) + 1
    offset = mod(((elements_in_head(T) - i) * bps) % UInt, 8 * sizeof(UInt))
    mask = UInt(1) << bps - 1
    return right_shift(@inbounds(seq.data[index]), offset) & mask
end

# This is usually type unstable, but in user code, users may use constant-folded ranges,
# e.g. f(x) = x[2:4]. In this case, we need it to compile to very efficient code.
# Hence, it MUST use @inline
@inline function Base.getindex(kmer::Kmer{A}, range::AbstractUnitRange{<:Integer}) where {A}
    @boundscheck checkbounds(kmer, range)
    K = length(range)
    iszero(K) && return Kmer{A, 0, 0}(unsafe, ())
    (i1, _) = BioSequences.bitindex(kmer, first(range))
    (i2, o2) = BioSequences.bitindex(kmer, last(range))
    data = kmer.data[i1:i2]
    (_, data) = rightshift_carry(data, o2, zero(UInt))
    T = derive_type(Kmer{A, K})
    N = nsize(T)
    # After the shift, the first coding element may be unused
    new_data = if N < length(data)
        tail(data)
    else
        data
    end
    return T(unsafe, (first(new_data) & get_mask(T), tail(new_data)...))
end

# Same as above: This needs to be able to inline if the indices are known statically
@inline function Base.getindex(kmer::Kmer{A}, indices::AbstractVector{Bool}) where {A}
    @boundscheck checkbounds(eachindex(kmer), indices)
    K = sum(indices)
    N = n_coding_elements(Kmer{A, K})
    T = Kmer{A, K, N}
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for (i, bool) in enumerate(indices)
        bool || continue
        (_, data) =
            leftshift_carry(data, nbits, BioSequences.extract_encoded_element(kmer, i))
    end
    return T(unsafe, data)
end

function Base.getindex(kmer::Kmer{A}, indices::AbstractVector{<:Integer}) where {A}
    K = length(indices)
    N = n_coding_elements(Kmer{A, K})
    T = Kmer{A, K, N}
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for i in indices
        checkbounds(kmer, i)
        (_, data) =
            leftshift_carry(data, nbits, BioSequences.extract_encoded_element(kmer, i))
    end
    return T(unsafe, data)
end

@inline function BioSequences.bitindex(kmer::Kmer, i::Integer)
    return BioSequences.bitindex(kmer, UInt(i)::UInt)
end

@inline function BioSequences.bitindex(kmer::Kmer, i::UInt)::Tuple{UInt, UInt}
    bps = BioSequences.bits_per_symbol(kmer) % UInt
    bpe = (8 * sizeof(UInt)) % UInt
    num = (UInt(i) - UInt(1) + n_unused(typeof(kmer)) % UInt) * bps
    (i, o) = divrem(num, bpe)
    o = bpe - o - bps
    return i + 1, o
end

@inline function Base.setindex(kmer::Kmer, i::Integer, s)
    @boundscheck checkbounds(kmer, i)
    bps = BioSequences.bits_per_symbol(kmer)
    iszero(bps) && return kmer
    symbol = convert(eltype(kmer), s)
    encoding = UInt(BioSequences.encode(Alphabet(kmer), symbol))
    (i, o) = BioSequences.bitindex(kmer, i % UInt)
    element = @inbounds kmer.data[i]
    mask = left_shift(UInt(1) << bps - 1, o)
    element &= ~mask
    element |= left_shift(encoding, o)
    return typeof(kmer)(unsafe, @inbounds Base.setindex(kmer.data, element, i))
end
