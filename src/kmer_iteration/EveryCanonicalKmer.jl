
struct EveryCanonicalKmer{T<:Kmer,S<:BioSequence{<:NucleicAcidAlphabet}} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
end

"""
    EveryCanonicalKmer(seq::BioSequence{A}, ::Val{K}) where {A<:NucleicAcidAlphabet,K}

Initialize an iterator over all overlapping k-mers in a sequence `seq` skipping
ambiguous nucleotides without changing the reading frame.
"""
function EveryCanonicalKmer(seq::BioSequence{A}, ::Val{K}) where {A<:NucleicAcidAlphabet,K}
    T′ = kmertype(Kmer{A,K})
    checkmer(T′) # Should inline and constant fold.
    return EveryCanonicalKmer{T′,typeof(seq)}(seq, 1, lastindex(seq))
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
        #return (pos, Kmer{A,K,N}(min(fwkmer, rvkmer))), (i, fwkmer, rvkmer)
    end
end