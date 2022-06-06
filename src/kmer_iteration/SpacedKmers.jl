

struct SpacedKmers{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

function SpacedKmers(::Type{T}, seq::S, step::Integer) where {T<:Kmer,S<:BioSequence}
    T′ = kmertype(T)
    checkmer(T′) # Should inline and constant fold.
    if step <= 1
        throw(ArgumentError("step size must be greater than 1"))
    end
    filled = max(0, ksize(T′) - step)
    increment = max(1, step - ksize(T′) + 1)
    return SpacedKmers{T′,typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
end

"""
    SpacedKmers(seq::BioSequence, ::Val{K}, step::Integer) where {K}

Initialize an iterator over k-mers separated by a `step` parameter, in a
sequence `seq` skipping ambiguous nucleotides without changing the reading frame.
"""
function SpacedKmers(seq::BioSequence{A}, ::Val{K}, step::Integer) where {A,K}
    T′ = kmertype(Kmer{A,K})
    checkmer(T′) # Should inline and constant fold.
    if step <= 1
        throw(ArgumentError("step size must be greater than 1"))
    end
    filled = max(0, K - step)
    increment = max(1, step - K + 1)
    return SpacedKmers{T′,typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
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
