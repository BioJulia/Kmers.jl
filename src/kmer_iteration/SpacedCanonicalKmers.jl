

struct SpacedCanonicalKmers{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

"""
    SpacedCanonicalKmers(seq::BioSequence, ::Val{K}, step::Integer) where {K}

Initialize an iterator over k-mers separated by a `step` parameter, in a
sequence `seq` skipping ambiguous nucleotides without changing the reading frame.
"""
function SpacedCanonicalKmers(seq::BioSequence{A}, ::Val{K}, step::Integer) where {A,K}
    T′ = kmertype(Kmer{A,K})
    checkmer(T′) # Should inline and constant fold.
    if step <= 1
        throw(ArgumentError("step size must be greater than 1"))
    end
    filled = max(0, K - step)
    increment = max(1, step - K + 1)
    return SpacedCanonicalKmers{T′,typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
end

Base.step(x::SpacedCanonicalKmers) = x.step

@inline function Base.iterate(it::SpacedCanonicalKmers{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    fwkmer = _build_kmer_data(Kmer{A,K,N}, it.seq, 1)
    if isnothing(fwkmer)
        return nothing
    else
        rshift = n_unused(Kmer{A,K,N}) * BioSequences.bits_per_symbol(A()) # Based on alphabet type, should constant fold.
        rvkmer = rightshift_carry(_reverse(BioSequences.BitsPerSymbol(A()), _complement_bitpar(A(), fwkmer...)...), rshift)
        return (1, min(Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer))), (K, fwkmer, rvkmer)
    end
end

@inline function Base.iterate(it::SpacedCanonicalKmers{Kmer{A,K,N},LongSequence{A}}, state) where {A,K,N}
    i, fwkmer, rvkmer = state
    filled = it.filled
    i += it.increment
    
    for _ in filled:K-1
        if i > it.stop
            return nothing
        else
            bps = BioSequences.bits_per_symbol(A()) # Based on type info, should constant fold.
            rshift = (64 - (n_unused(Kmer{A,K,N}) + 1) * bps) # Based on type info, should constant fold.
            mask = (one(UInt64) << bps) - one(UInt64) # Based on type info, should constant fold.
            fbits = UInt64(BioSequences.extract_encoded_element(it.seq, i))
            rbits = (BioSequences.complement_bitpar(fbits, A()) & mask) << rshift
            fwkmer = leftshift_carry(fwkmer, bps, fbits)
            rvkmer = rightshift_carry(rvkmer, bps, rbits)
            i += 1
        end
    end
    pos = i - K + 1
    return (pos, min(Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer))), (i, fwkmer, rvkmer)
end
#=
@inline function Base.iterate(it::SpacedCanonicalKmers{Kmer{A,K,N},LongSequence{A}}, state = (it.start - it.increment, 1, 0, blank_ntuple(Kmer{A,K,N}), blank_ntuple(Kmer{A,K,N}))
    ) where {A,K,N}
    i, pos, filled, fwkmer, rvkmer = state
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(UInt8, @inbounds getindex(it.seq, i))
        @inbounds fbits = UInt64(kmerbits[nt + 1])
        rbits = ~fbits & typeof(fbits)(0x03)
        if fbits == 0xff # ambiguous
            filled = 0
            # Find the beginning of next possible kmer after i
            pos = i + it.step - Core.Intrinsics.urem_int(i - pos, it.step)
            i = pos - 1
        else
            filled += 1
            fwkmer = leftshift_carry(fwkmer, 2, fbits)
            rvkmer = rightshift_carry(rvkmer, 2, UInt64(rbits) << (62 - (64N - 2K)))
        end
        if filled == K
            state = (i, i - K + 1 + it.step, it.filled, fwkmer, rvkmer)
            return KmerAt(pos, Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer)), state
        end
        i += 1
    end
    return nothing
end
=#