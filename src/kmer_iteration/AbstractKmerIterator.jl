###
### Kmer Iteration
###
### Abstract Kmer Iterator type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

### Type for storing the result of Kmer iteration.

"""
A type yielded by AbstractKmerIterators

Represents the kmer at a given position - both on the forward and reverse strand.
"""
struct KmerAt{T<:Kmer}
    position::Int
    fw::T
    bw::T
end

@inline function Base.iterate(x::KmerAt, state = 1)
    if 1 ≤ state ≤ 3
        return getfield(x, state), state + 1
    else
        return nothing
    end
end

@inline position(x::KmerAt) = x.position
@inline fwmer(x::KmerAt) = x.fw
@inline bwmer(x::KmerAt) = x.bw
@inline canonical(x::KmerAt) = min(fwmer(x), bwmer(x))

function Base.show(io::IO, ::MIME"text/plain", x::KmerAt)
    println(io, "Kmer at position:", x.position)
    print(io, "Forward: ")
    showcompact(io, fwmer(x))
    print(io, "\nBackward: ")
    showcompact(io, bwmer(x))
    print(io, '\n')
end

abstract type AbstractKmerIterator{T<:Kmer,S<:BioSequence} end
#abstract type CanonicalKmerIterator{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}  end

@inline Base.eltype(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Tuple{UInt64,T}
#@inline Base.eltype(::Type{<:CanonicalKmerIterator{T,S}}) where {T,S} = KmerAt{T}

@inline Base.IteratorSize(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Base.HasLength()

@inline function Base.length(it::AbstractKmerIterator{T,S}) where {T,S}
    return max(0, fld(it.stop - it.start + 1 - ksize(T), step(it)) + 1)
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