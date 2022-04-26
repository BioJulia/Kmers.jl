###
### Kmer Iteration
###
### Iterator type over every kmer in a sequence - overlapping.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


struct EveryKmer{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
end

"""
    everykmer(::Val{K}, seq::BioSequence) where {K}

Initialize an iterator over all overlapping k-mers in a sequence `seq` skipping
ambiguous nucleotides without changing the reading frame.
"""
function EveryKmer(seq::BioSequence{A}, ::Val{K}) where {A,K}
    T′ = kmertype(Kmer{A,K})
    checkmer(T′) # Should inline and constant fold.
    return EveryKmer{T′,typeof(seq)}(seq, 1, lastindex(seq))
end

Base.step(x::EveryKmer) = 1

## Initial iteration without state.
@inline function Base.iterate(it::EveryKmer{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    kmer = _build_kmer_data(Kmer{A,K,N}, it.seq, 1)
    if isnothing(kmer)
        return nothing
    else
        # Get the reverse.
        alph = Alphabet(Kmer{A,K,N})
        return (1, Kmer{A,K,N}(kmer)), (K, kmer)
    end
end

@inline function Base.iterate(it::EveryKmer{Kmer{A,K,N},LongSequence{A}}, state) where {A,K,N}
    i, fwkmer = state
    i += 1
    if i > it.stop
        return nothing
    else
        bps = BioSequences.bits_per_symbol(A()) # Based on type info, should constant fold.
        bits = UInt64(BioSequences.extract_encoded_element(it.seq, i))
        kmer = leftshift_carry(fwkmer, bps, bits)
        pos = i - K + 1
        return (pos, Kmer{A,K,N}(kmer)), (i, kmer)
    end
end