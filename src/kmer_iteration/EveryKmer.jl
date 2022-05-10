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

function EveryKmer(::Type{T}, seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
    return EveryKmer{T,S}(seq, start, stop)
end

"""
    EveryKmer(::Val{K}, seq::BioSequence) where {K}

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

## Special case where iterating over 2-Bit encoded kmers in a 4-Bit encoded sequence,
## behaviour is to produce kmers by skipping over the ambiguous sites.

const kmerbits = (UInt64(0xff), UInt64(0x00), UInt64(0x01), UInt64(0xff),
                  UInt64(0x02), UInt64(0xff), UInt64(0xff), UInt64(0xff),
                  UInt64(0x03), UInt64(0xff), UInt64(0xff), UInt64(0xff),
                  UInt64(0xff), UInt64(0xff), UInt64(0xff), UInt64(0xff))

@inline Base.IteratorSize(::Type{<:EveryKmer{Kmer{A,N,K},S}}) where {A<:NucleicAcidAlphabet{2},N,K,B<:NucleicAcidAlphabet{4},S<:BioSequence{B}} = Base.SizeUnknown()

@inline function Base.iterate(it::EveryKmer{Kmer{A,K,N},S},
                              state = (it.start - 1, 1, blank_ntuple(Kmer{A,K,N}))
                              ) where {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4},S<:BioSequence{B},K,N}
    
    i, filled, fwkmer = state
    i += 1
    filled -= 1

    while i ≤ it.stop
        @inbounds nt = reinterpret(UInt8, it.seq[i])
        @inbounds fbits = kmerbits[nt + 1]
        fwkmer = leftshift_carry(fwkmer, 2, fbits)
        filled = ifelse(fbits == UInt64(0xff), 0, filled + 1)
        if filled == K
            return (i - K + 1, Kmer{A,K,N}(fwkmer)), (i, filled, fwkmer)
        end
        i += 1
    end
    return nothing
end