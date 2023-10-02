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

    function SpacedKmers{A, K, St, S}(seq) where {A, K, St, S}
        (K isa Int && K > 0) || error("K must be an Int > 0")
        (St isa Int && St > 0) || error("St must be an Int > 0")
        new{A, K, St, S}(seq)
    end
end

source_type(::Type{SpacedKmers{A, K, St, S}}) where {A, K, St, S} = S

function Base.length(it::SpacedKmers{A, K, St}) where {A, K, St}
    Base.IteratorSize(typeof(it)) == Base.HasLength() || throw(MethodError(length, (it,)))
    available_starting_positions = length(it.seq) - ksize(eltype(it)) + 1
    div(available_starting_positions, St)
end

# Constructors
SpacedKmers{A, K, St}(seq) where {A, K, St} = SpacedKmers{A, K, St, typeof(seq)}(seq)

# Iterators
function Base.iterate(it::SpacedKmers{A, K, St, <:BioSequence{A}}) where {A, K, St}
    if St ≥ K
        iterate_copy_nomask(it, 1)
    else
        x = iterate_copy_nomask(it, 1)
        x === nothing && return nothing
        (kmer, _) = x
        return (kmer, (K + 1, kmer.data))
    end
end

function Base.iterate(it::SpacedKmers{A, K, St, <:BioSequence{A}}, state) where {A, K, St}
    if St ≥ K
        iterate_copy_nomask(it, state)
    else
        iterate_copy_mask(it, state)
    end
end

# Called when St ≥ K, and the encoding in seq matches that of the kmer.
# We can build the kmer from scratch at every iteration, simplifying the code
@inline function iterate_copy_nomask(it::SpacedKmers{A, K, St}, state::Int) where {A, K, St}
    seq = it.seq
    len = length(seq)
    bps = BioSequences.bits_per_symbol(A())
    remaining = K
    data = zero_tuple(eltype(it))
    while true
        state > len && return nothing
        encoding = UInt(BioSequences.extract_encoded_element(seq, state))
        (_, data) = leftshift_carry(data, bps, encoding)
        state += 1
        remaining -= 1
        iszero(remaining) && return (eltype(it)(unsafe, data), state + max(0, St - K))
    end
end

# Called when St < K, and the encoding in seq matches that of the kmer.
# We can copy the encoding right over, and we need to preserve some data in the kmer
# between iterations
@inline function iterate_copy_mask(
    it::SpacedKmers{A, K, St},
    state::Tuple{Int, Tuple{Vararg{UInt}}},
) where {A, K, St}
    seq = it.seq
    len = length(seq)
    bps = BioSequences.bits_per_symbol(A())
    remaining = St
    (index, data) = state
    while true
        index > len && return nothing
        encoding = UInt(BioSequences.extract_encoded_element(seq, index))
        (_, data) = leftshift_carry(data, bps, encoding)
        index += 1
        remaining -= 1
        if iszero(remaining)
            # Mask out unused bits before we return the kmer.
            (head, rest...) = data
            kmer = eltype(it)(unsafe, (head & get_mask(eltype(it)), rest...))
            return (kmer, (index, data))
        end
    end
end

# TODO: Methods:
# 4 -> 2 bit
# 2 -> 4 bit?
# Byte sequence: 2 bit
# Byte sequence: other alphabets
