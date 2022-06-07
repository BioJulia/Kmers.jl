###
### Kmer Iteration
###
### Iterator type over every kmer in a sequence - overlapping.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
An iterator over every valid overlapping `T<:Kmer` in a given longer
`BioSequence` between a `start` and `stop` position. 

!!! note
    Typically, the alphabet of the Kmer type matches the alphabet of the input
    BioSequence. In these cases, the iterator will have `Base.IteratorSize` of
    `Base.HasLength`, and successive kmers produced by the iterator will overlap
    by K - 1 bases.
    
    However, in the specific case of iterating over kmers in a DNA or RNA sequence, you
    may iterate over a Kmers where the alphabet is a NucleicAcidAlphabet{2}, but
    the input BioSequence has a NucleicAcidAlphabet{4}.
    
    In this case then the iterator will skip over positions in the BioSequence
    with characters that are not supported by the Kmer type's NucleicAcidAlphabet{2}.
    
    As a result, the overlap between successive kmers may not reliably be K - 1,
    and the iterator will have `Base.IteratorSize` of `Base.SizeUnknown`.
"""
struct EveryKmer{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
    
    function EveryKmer{T,S}(seq::S, start::Int = firstindex(seq), stop::Int = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
        T′ = kmertype(T)
        checkmer(T′) # Should inline and constant fold.
        return new{T′,S}(seq, start, stop)
    end
end

"""
    EveryKmer{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}

Convenience outer constructor so you don't have to specify `S` along with `T`.

E.g. Instead of `EveryKmer{DNACodon,typeof(s)}(s)`, you can just use `EveryKmer{DNACodon}(s)`
"""
function EveryKmer{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
    return EveryKmer{T,S}(seq, start, stop)
end

"""
    EveryKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}

Convenience outer constructor so yyou don't have to specify full `Kmer` typing.

In order to deduce `Kmer{A,K,N}`, `A` is taken from the input `seq` type, `K` is
taken from `::Val{K}`, and `N` is deduced using `A` and `K`.

E.g. Instead of `EveryKmer{DNAKmer{3,1}}(s)`, or `EveryKmer{DNACodon}(s)`,
you can use `EveryKmer(s, Val(3))`
"""
function EveryKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
    return EveryKmer{Kmer{A,K}}(seq, start, stop)
end

Base.step(x::EveryKmer) = 1

## Initial iteration without state.
@inline function Base.iterate(it::EveryKmer{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    kmer = _build_kmer_data(Kmer{A,K,N}, it.seq, 1)
    if isnothing(kmer)
        return nothing
    else
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