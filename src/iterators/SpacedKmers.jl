
"""
    SpacedKmers{T,S}(seq::S, step::Int, start::Int, stop::Int) where {T<:Kmer,S<:BioSequence}

An iterator over every valid `T<:Kmer` separated by a `step` parameter, in a given
longer `BioSequence`, between a `start` and `stop` position.

!!! note
    Typically, the alphabet of the Kmer type matches the alphabet of the input
    BioSequence. In these cases, the iterator will have `Base.IteratorSize` of
    `Base.HasLength`, and successive kmers produced by the iterator will overlap
    by `max(0, K - step)` bases.
    
    However, in the specific case of iterating over kmers in a DNA or RNA sequence, you
    may iterate over a Kmers where the alphabet is a NucleicAcidAlphabet{2}, but
    the input BioSequence has a NucleicAcidAlphabet{4}.
    
    In this case then the iterator will skip over positions in the BioSequence
    with characters that are not supported by the Kmer type's NucleicAcidAlphabet{2}.
    
    As a result, the overlap between successive kmers may not consistent, but the
    reading frame will be preserved.
    In addition, the iterator will have `Base.IteratorSize` of `Base.SizeUnknown`.
"""
struct SpacedKmers{A <: Alphabet, K, St, S} <: AbstractKmerIterator{A, K}
    seq::S

    function SpacedKmer{A, K, St, S}(seq) where {A, K, St, S}
        K isa Int || K > 0 || error("K must be an Int > 0")
        St isa Int || St > 0 ||  error("St must be an Int > 0")
        new{A, K, St, S}(seq)
    end
end

# Constructors

# Iterators
function Base.iterate(
    it::SpacedKmers{A, K, St, <:BioSequence{A}},
    state=(1, zero_tuple(eltype(it)), K)
) where {A, K, St}
    iterate_copy(it, state)
end

# Just copy the encoding straight over
@inline function iterate_copy(it::SpacedKmers{A, K, St}, state::Tuple{Int, <:Tuple{Vararg{UInt}}, Int}) where {A, K, St}
    (index, data, remaining) = state
    original_index = index
    seq = it.seq
    len = length(seq)
    bps = BioSequences.bits_per_symbol(A())
    # TODO: Don't double check remaining and index
    while !iszero(remaining) && index ≤ len
        encoding = UInt(BioSequences.extract_encoded_element(seq, index))
        (_, data) = leftshift_carry(data, bps, encoding)
    end
    return if iszero(remaining)
        kmer = eltype(it)(unsafe, data)
        remaining = min(K, St)
        state = (
            original_index + St,
            # data, left shift carry clear elements
            remaining
        )
        (kmer, state)
    else
        nothing
    end
end