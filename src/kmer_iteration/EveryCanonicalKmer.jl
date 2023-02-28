"""
    EveryCanonicalKmer{T,S}(seq::S, start::Int = firstindex(seq), stop::Int = lastindex(seq)) where {T<:Kmer,S<:BioSequence}

An iterator over every canonical valid overlapping `T<:Kmer` in a given longer 
`BioSequence`, between a `start` and `stop` position.

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
struct EveryCanonicalKmer{T<:Kmer,S<:BioSequence{<:NucleicAcidAlphabet}} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
    
    function EveryCanonicalKmer{T,S}(seq::S, start::Int = firstindex(seq), stop::Int = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
        T′ = kmertype(T)
        checkmer(T′) # Should inline and constant fold.
        return new{T′,S}(seq, start, stop)
    end
end

"""
    EveryCanonicalKmer{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}

Convenience outer constructor so you don't have to specify `S` along with `T`.

E.g. Instead of `EveryCanonicalKmer{DNACodon,typeof(s)}(s)`, you can just use `EveryCanonicalKmer{DNACodon}(s)`
"""
function EveryCanonicalKmer{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
    return EveryCanonicalKmer{T,S}(seq, start, stop)
end

"""
    EveryCanonicalKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}

Convenience outer constructor so yyou don't have to specify full `Kmer` typing.

In order to deduce `Kmer{A,K,N}`, `A` is taken from the input `seq` type, `K` is
taken from `::Val{K}`, and `N` is deduced using `A` and `K`.

E.g. Instead of `EveryCanonicalKmer{DNAKmer{3,1}}(s)`, or `EveryCanonicalKmer{DNACodon}(s)`,
you can use `EveryCanonicalKmer(s, Val(3))`
"""
function EveryCanonicalKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
    return EveryCanonicalKmer{Kmer{A,K}}(seq, start, stop)
end

Base.step(x::EveryCanonicalKmer) = 1



## Initial iteration without state.
@inline function Base.iterate(it::EveryCanonicalKmer{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    fwkmer = _build_kmer_data(Kmer{A,K,N}, it.seq, it.start)
    if isnothing(fwkmer)
        return nothing
    else
        rshift = n_unused(Kmer{A,K,N}) * BioSequences.bits_per_symbol(A()) # Based on alphabet type, should constant fold.
        rvkmer = rightshift_carry(_reverse(BioSequences.BitsPerSymbol(A()), _complement_bitpar(A(), fwkmer...)...), rshift)
        return (it.start, Kmer{A,K,N}(min(fwkmer, rvkmer))), (it.start + K - 1, fwkmer, rvkmer)
    end
end

@inline function Base.iterate(it::EveryCanonicalKmer{Kmer{A,K,N},LongSequence{A}}, state) where {A,K,N}
    i, fwkmer, rvkmer = state
    i += 1
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
        pos = i - K + 1
        return (pos, min(Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer))), (i, fwkmer, rvkmer)
    end
end

@inline Base.IteratorSize(::Type{<:EveryCanonicalKmer{Kmer{A,N,K},LongSequence{B}}}) where {A<:NucleicAcidAlphabet{2},N,K,B<:NucleicAcidAlphabet{4}} = Base.SizeUnknown()

@inline function Base.iterate(it::EveryCanonicalKmer{Kmer{A,K,N},LongSequence{B}},
                              state = (it.start - 1, 1, blank_ntuple(Kmer{A,K,N}), blank_ntuple(Kmer{A,K,N}))
                              ) where {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4},K,N}
    
    i, filled, fwkmer, rvkmer = state
    i += 1
    filled -= 1

    rshift = (64 - (n_unused(Kmer{A,K,N}) + 1) * 2) # Based on type info, should constant fold.
    mask = (one(UInt64) << 2) - one(UInt64) # Based on type info, should constant fold.

    while i ≤ it.stop
        @inbounds nt = reinterpret(UInt8, it.seq[i])
        @inbounds fbits = kmerbits[nt + 1]
        rbits = (BioSequences.complement_bitpar(fbits, A()) & mask) << rshift
        fwkmer = leftshift_carry(fwkmer, 2, fbits)
        rvkmer = rightshift_carry(rvkmer, 2, rbits)
        filled = ifelse(fbits == UInt64(0xff), 0, filled + 1)
        if filled == K
            return (i - K + 1, min(Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer))), (i, filled, fwkmer, rvkmer)
        end
        i += 1
    end
    return nothing
end