@inline BioSequences.encoded_data_eltype(::Type{<:Kmer}) = UInt64

@inline function BioSequences.extract_encoded_element(seq::Kmer, i::Integer)
    bi = BioSequences.bitindex(seq, i % UInt)
    return BioSequences.extract_encoded_element(bi, seq.data)
end

@inline Base.copy(seq::Kmer) = typeof(seq)(seq.data)

@inline encoded_data(x::Kmer) = x.data

@inline BioSequences.bitindex(seq::Kmer, i::Integer) = BioSequences.bitindex(BioSequences.BitsPerSymbol(seq), BioSequences.encoded_data_eltype(typeof(seq)), i + n_unused(seq))


"""
Base.getindex(seq::Kmer, i::UnitRange)

Slice a Kmer by a UnitRange.

!!! warning
    Using this function will introduce performance penalties in your code if
    you pass values of `i` that are not constants that can be propagated.
"""
@inline function Base.getindex(seq::Kmer{A}, i::UnitRange) where A
    @boundscheck Base.checkbounds(seq, i)
    ind(s, i) = BioSequences.index(BioSequences.bitindex(s, i))
    off(s, i) = BioSequences.offset(BioSequences.bitindex(s, i))
    isempty(i) && return Kmer{A, 0, 0}(())
    rshift = (64 - off(seq, last(i) + 1)) & 63
    stop = ind(seq, last(i))
    start = BioSequences.index(BioSequences.bitindex(seq, first(i)) + rshift)
    data = Kmers.rightshift_carry(seq.data, rshift)
    T = Kmers.kmertype(Kmer{A, length(i)})
    return T(data[start:stop])
end