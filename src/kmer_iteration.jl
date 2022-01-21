###
### Kmer Iteration
###
### Iterator over all k-mers in a biological sequence.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Note about the variable names:

# it.step is the distance from start of one kmer to start of next
#
# filled is the number of nucleotides in a kmer that has the correct value set.
# e.g. when moving from xxxxxxx to yyyyyyy, filled goes from K to 3.
# when moving from xxxxxxx to zzzzzzz, filled goes from K to 0.
#
# increment is how far the index jumps ahead when going to the next kmer.
# for close kmers where no jump is possible, it's 1, else it can be arbitrary high.
#
# For a kmeriterator that first emits xxxxxxx, then zzzzzzz:
#
#       |--------- step ---------|
#       xxxxxxx                  zzzzzzz
# -------------------------------------------------------------
#           yyyyyyy
#             |---- increment ---|
#
# The state returned at each iteration is the state upon return, not the state
# needed for the following iteration.

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

struct EveryKmerIterator{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
end

struct SpacedKmerIterator{T<:Kmer,S<:BioSequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

"""
    kmers(::Type{T}, seq::BioSequence) where {T<:Kmer}

Initialize an iterator over all overlapping k-mers in a sequence `seq` skipping
ambiguous nucleotides without changing the reading frame.
"""
function kmers(::Type{T}, seq::BioSequence) where {T<:Kmer}
    T′ = kmertype(T)
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    end
    checkmer(T′) # Should inline and constant fold.
    return EveryKmerIterator{T′,typeof(seq)}(seq, 1, lastindex(seq))
end

"""
    kmers(::Type{T}, seq::BioSequence, step::Integer) where {T<:Kmer}

Initialize an iterator over k-mers separated by a `step` parameter, in a
sequence `seq` skipping ambiguous nucleotides without changing the reading frame.
"""
function kmers(::Type{T}, seq::BioSequence, step::Integer) where {T<:Kmer}
    T′ = kmertype(T)
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    checkmer(T′) # Should inline and constant fold.
    if step == 1
        return EveryKmerIterator{T′,typeof(seq)}(seq, 1, lastindex(seq))
    else
        filled = max(0, ksize(T′) - step)
        increment = max(1, step - ksize(T′) + 1)
        return SpacedKmerIterator{T′,typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
    end
end

@inline Base.eltype(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = KmerAt{T}

@inline Base.IteratorSize(::Type{<:AbstractKmerIterator{<:Kmer,LongSequence{<:NucleicAcidAlphabet{4}}}}) = Base.SizeUnknown()
@inline Base.IteratorSize(::Type{<:AbstractKmerIterator{<:Kmer,LongSequence{<:NucleicAcidAlphabet{2}}}}) = Base.HasLength()

@inline function Base.length(it::AbstractKmerIterator{T,LongSequence{<:NucleicAcidAlphabet{2}}}) where {T<:Kmer}
    return max(0, fld(it.stop - it.start + 1 - ksize(T), step(it)) + 1)
end

Base.step(x::EveryKmerIterator) = 1
Base.step(x::SpacedKmerIterator) = x.step

# Initial iteration where the alphabet between seq and kmer match.
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

@inline function Base.iterate(it::EveryKmerIterator{Kmer{A,K,N},LongSequence{A}}, state) where {A,K,N}
    i, fwkmer, rvkmer = state
    i += 1
    if i > it.stop
        return nothing
    else
        bps = BioSequences.bits_per_symbol(A()) # Based on type info, should constant fold.
        rshift = (62 - n_unused(Kmer{A,K,N}) * bps) # Based on type info, should constant fold.
        mask = (one(UInt64) << bps) - one(UInt64) # Based on type info, should constant fold.
        
        fbits = UInt64(BioSequences.extract_encoded_element(it.seq, i))
        rbits = (BioSequences.complement_bitpar(fbits, A()) & mask) << rshift
        fwkmer = leftshift_carry(fwkmer, bps, fbits)
        rvkmer = rightshift_carry(rvkmer, bps, rbits)
        pos = i - K + 1
        return KmerAt(pos, Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer)), (i, fwkmer, rvkmer)
    end
end

@inline function Base.iterate(
    it::EveryKmerIterator{Kmer{A,K,N},LongSequence{<:NucleicAcidAlphabet{4}}},
    state = (it.start - 1, 1, blank_ntuple(Kmer{A,K,N}), blank_ntuple(Kmer{A,K,N}))
    ) where {A<:NucleicAcidAlphabet{2},K,N}
    
    i, filled, fkmer, rkmer = state
    i += 1
    filled -= 1
    rshift = (62 - (64N - 2K))
    
    while i ≤ it.stop
        nt = reinterpret(UInt8, @inbounds it.seq[i])
        # Do this rather than `encode` or decode to avoid throwing an exception,
        # and instead just skip over those bases.
        fbits = UInt64(@inbounds twobitnucs[nt + 0x01])
        rbits = (~fbits & UInt64(0x03)) << rshift
        fkmer = leftshift_carry(fkmer, 2, fbits)
        rkmer = rightshift_carry(rkmer, 2, )
        filled = ifelse(fbits == 0xff, 0, filled + 1)
        if filled == K
            pos = i - K + 1
            return KmerAt(pos, Kmer{A,K,N}(fkmer), Kmer{A,K,N}(rkmer)), (i, K, fkmer, rkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedKmerIterator{Kmer{A,K,N},S}, state=(it.start - it.increment, 1, 0, blank_ntuple(Kmer{A,K,N}), blank_ntuple(Kmer{A,K,N}))
    ) where {A,K,N,S<:LongSequence{<:NucleicAcidAlphabet{4}}}
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
