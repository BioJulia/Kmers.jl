
"""
An iterator over every valid `T<:Kmer` separated by a `step` parameter, in a given
longer `BioSequence`.

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
struct SpacedKmers{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
    
    function SpacedKmers{T,S}(seq::S, step::Int, start::Int, stop::Int) where {T<:Kmer,S<:BioSequence}
        T′ = kmertype(T)
        checkmer(T′) # Should inline and constant fold.
        if step <= 1
            throw(ArgumentError("step size must be greater than 1"))
        end
        filled = max(0, ksize(T′) - step)
        increment = max(1, step - ksize(T′) + 1)
        return new{T′,S}(seq, start, step, stop, filled, increment)
    end
end

"""
    SpacedKmers{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}

Convenience outer constructor so you don't have to specify `S` along with `T`.

E.g. Instead of `SpacedKmers{DNACodon,typeof(s)}(s, 3)`, you can just use `SpacedKmers{DNACodon}(s, 3)`
"""
function SpacedKmers{T}(seq::S, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
    return SpacedKmers{T,S}(seq, step, start, stop)
end

"""
    SpacedKmers(seq::BioSequence{A}, ::Val{K}, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {A,K}

Convenience outer constructor so yyou don't have to specify full `Kmer` typing.

In order to deduce `Kmer{A,K,N}`, `A` is taken from the input `seq` type, `K` is
taken from `::Val{K}`, and `N` is deduced using `A` and `K`.

E.g. Instead of `SpacedKmers{DNAKmer{3,1}}(s, 3)`, or `SpacedKmers{DNACodon}(s, 3)`,
you can use `SpacedKmers(s, Val(3), 3)`
"""
function SpacedKmers(seq::BioSequence{A}, ::Val{K}, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
    return SpacedKmers{Kmer{A,K}}(seq, step, start, stop)
end

Base.step(x::SpacedKmers) = x.step

@inline function Base.iterate(it::SpacedKmers{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    kmer = _build_kmer_data(Kmer{A,K,N}, it.seq, 1)
    if isnothing(kmer)
        return nothing
    else
        # Get the reverse.
        alph = Alphabet(Kmer{A,K,N})
        return (1, Kmer{A,K,N}(kmer)), (K, kmer)
    end
end

@inline function Base.iterate(it::SpacedKmers{Kmer{A,K,N},LongSequence{A}}, state) where {A,K,N}
    i, kmer = state
    filled = it.filled
    i += it.increment
    
    for _ in filled:K-1
        if i > it.stop
            return nothing
        else
            bps = BioSequences.bits_per_symbol(A()) # Based on type info, should constant fold.
            bits = UInt64(BioSequences.extract_encoded_element(it.seq, i))
            kmer = leftshift_carry(kmer, bps, bits)
            i += 1
        end
    end
    pos = i - K + 1
    return (pos, Kmer{A,K,N}(kmer)), (i, kmer)
end

@inline function Base.iterate(it::SpacedKmers{Kmer{A,K,N},LongSequence{B}}, state = (it.start - it.increment, 1, 0, blank_ntuple(Kmer{A,K,N}))
    ) where {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4},K,N}
    i, pos, filled, kmer = state
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(UInt8, @inbounds getindex(it.seq, i))
        @inbounds bits = UInt64(kmerbits[nt + 1])
        if bits == 0xff # ambiguous
            filled = 0
            # Find the beginning of next possible kmer after i
            pos = i + it.step - Core.Intrinsics.urem_int(i - pos, it.step)
            i = pos - 1
        else
            filled += 1
            kmer = leftshift_carry(kmer, 2, bits)
        end
        if filled == K
            state = (i, i - K + 1 + it.step, it.filled, kmer)
            return (pos, Kmer{A,K,N}(kmer)), state
        end
        i += 1
    end
    return nothing
end
