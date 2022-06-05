###
### Kmer Iteration
###
### Abstract Kmer Iterator type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

### Type for storing the result of Kmer iteration.

abstract type AbstractKmerIterator{T<:Kmer,S<:BioSequence} end

@inline Base.eltype(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Tuple{UInt64,T}

@inline Base.IteratorSize(::Type{<:AbstractKmerIterator{Kmer{A,K,N},S}}) where {A,S<:BioSequence{A},K,N} = Base.HasLength()
@inline Base.IteratorSize(::Type{<:AbstractKmerIterator{Kmer{A,K,N},S}}) where {A,B,S<:BioSequence{B},K,N} = Base.SizeUnknown()

@inline function Base.length(it::AbstractKmerIterator{Kmer{A,K,N},S}) where {A,K,N,S<:BioSequence{A}}
    return max(0, fld(it.stop - it.start + 1 - K, step(it)) + 1)
end

# Iteration where the Kmer and Seq alphabets match:

## Initial iteration without state.
@inline function Base.iterate(it::AbstractKmerIterator{Kmer{A,K,N},LongSequence{A}}) where {A,K,N}
    fwkmer = _build_kmer_data(Kmer{A,K,N}, it.seq, 1)
    if isnothing(fwkmer)
        return nothing
    else
        # Get the reverse.
        alph = Alphabet(Kmer{A,K,N})
        rshift = n_unused(Kmer{A,K,N}) * BioSequences.bits_per_symbol(alph) # Based on alphabet type, should constant fold.
        rvkmer = rightshift_carry(_reverse(BioSequences.BitsPerSymbol(alph), _complement_bitpar(alph, fwkmer...)...), rshift)
        return KmerAt{Kmer{A,K,N}}(1, Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer)), (K, fwkmer, rvkmer)
    end
end