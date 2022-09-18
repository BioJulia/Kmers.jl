###
### Kmer Type definition
###

# Include some basic tuple bitflipping ops - the secret sauce to efficiently
# manipping Kmer's static data. 
include("tuple_bitflipping.jl")

"""
    Kmers.Unsafe

Trait object used to access unsafe methods of functions.
`unsafe` is the singleton of `Unsafe`.
"""
struct Unsafe end
const unsafe = Unsafe()

"""
    Kmer{A<:Alphabet,K,N} <: BioSequence{A}

A parametric, immutable, bitstype for representing Kmers - short sequences.
Given the number of Kmers generated from raw sequencing reads, avoiding
repetetive memory allocation and triggering of garbage collection is important,
as is the ability to effectively pack Kmers into arrays and similar collections.

In practice that means we an immutable bitstype as the internal representation
of these sequences. Thankfully, this is not much of a limitation - kmers are
rarely manipulated and so by and large don't have to be mutable.

Excepting their immutability, they fulfill the rest of the API and behaviours
expected from a concrete `BioSequence` type, and non-mutating transformations
of the type are still defined.

!!! warning
    Given their immutability, `setindex` and mutating sequence transformations
    are not implemented for Kmers e.g. `reverse_complement!`. 
!!! tip
    Note that some sequence transformations that are not mutating are
    available, since they can return a new kmer value as a result e.g.
    `reverse_complement`. 
"""
struct Kmer{A<:Alphabet,K,N} <: BioSequence{A}
    data::NTuple{N,UInt64}
    
    # This unsafe method do not clip the head
    Kmer{A,K,N}(::Unsafe, data::NTuple{N,UInt64}) where {A<:Alphabet,K,N} = new{A,K,N}(data)

    function Kmer{A,K,N}(data::NTuple{N,UInt64}) where {A<:Alphabet,K,N}
        checkmer(Kmer{A,K,N})
        x = n_unused(Kmer{A,K,N}) * BioSequences.bits_per_symbol(A()) 
        return new(_cliphead(x, data...))
    end
end

BioSequences.encoded_data(seq::Kmer{A,K,N}) where {A,K,N} = seq.data

# Create a blank ntuple of appropriate length for a given Kmer with N.
@inline blank_ntuple(::Type{Kmer{A,K,N}}) where {A,K,N} = ntuple(x -> zero(UInt64), Val{N}())

# Aliases
"Shortcut for the type `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K,N} = Kmer{DNAAlphabet{2},K,N}

"Shortcut for the type `DNAKmer{27,1}`"
const DNA27mer = DNAKmer{27,1}

"Shortcut for the type `DNAKmer{31,1}`"
const DNA31mer = DNAKmer{31,1}

"Shortcut for the type `DNAKmer{63,2}`"
const DNA63mer = DNAKmer{63,2}

"Shortcut for the type `Kmer{RNAAlphabet{2},K,N}`"
const RNAKmer{K,N} = Kmer{RNAAlphabet{2},K,N}

"Shortcut for the type `RNAKmer{27,1}`"
const RNA27mer = RNAKmer{27,1}

"Shortcut for the type `RNAKmer{31,1}`"
const RNA31mer = RNAKmer{31,1}

"Shortcut for the type `RNAKmer{63,2}`"
const RNA63mer = RNAKmer{63,2}

"Shortcut for the type `Kmer{AminoAcidAlphabet,K,N}`"
const AAKmer{K,N} = Kmer{AminoAcidAlphabet,K,N}

"Shorthand for `DNAKmer{3,1}`"
const DNACodon = DNAKmer{3,1}

"Shorthand for `RNAKmer{3,1}`"
const RNACodon = RNAKmer{3,1}





include("constructors.jl")




"""
    kmertype(::Type{Kmer{A,K}}) where {A,K}
Resolve and incomplete kmer typing, computing the N parameter of
`Kmer{A,K,N}`, given only `Kmer{A,K}`.
## Example
```julia
julia> DNAKmer{63}
Kmer{DNAAlphabet{2},63,N} where N
julia> kmertype(DNAKmer{63})
Kmer{DNAAlphabet{2},63,2}
```
"""
@inline function kmertype(::Type{Kmer{A,K}}) where {A,K}
    return Kmer{A,K,BioSequences.seq_data_len(A, K)}
end
@inline kmertype(::Type{Kmer{A,K,N}}) where {A,K,N} = Kmer{A,K,N}

###
### Base Functions
###

@inline ksize(::Type{Kmer{A,K,N}}) where {A,K,N} = K
@inline nsize(::Type{Kmer{A,K,N}}) where {A,K,N} = N
@inline per_word_capacity(::Type{Kmer{A,K,N}}) where {A,K,N} = div(64, BioSequences.bits_per_symbol(A()))
@inline per_word_capacity(seq::Kmer) = per_word_capacity(typeof(seq))
@inline capacity(::Type{Kmer{A,K,N}}) where {A,K,N} = per_word_capacity(Kmer{A,K,N}) * N
@inline capacity(seq::Kmer) = capacity(typeof(seq))
@inline n_unused(::Type{Kmer{A,K,N}}) where {A,K,N} = capacity(Kmer{A,K,N}) - K
@inline n_unused(seq::Kmer) = n_unused(typeof(seq))
@inline elements_in_head(::Type{Kmer{A,K,N}}) where {A,K,N} = per_word_capacity(Kmer{A,K,N}) - n_unused(Kmer{A,K,N})
@inline elements_in_head(seq::Kmer) = elements_in_head(typeof(seq))

"""
    checkmer(::Type{Kmer{A,K,N}}) where {A,K,N}

Internal method - enforces good kmer type parameterisation.

For a given Kmer{A,K,N} of length K, the number of words used to
represent it (N) should be the minimum needed to contain all K symbols,
no larger (wasteful) no smaller (just... wrong).

Because it is used on type parameters / variables, these conditions should be
checked at compile time, and the branches / error throws eliminated when the
parameterisation of the Kmer type is good. 
"""
@inline function checkmer(::Type{Kmer{A,K,N}}) where {A,K,N}
    if K < 1
        throw(ArgumentError("Bad kmer parameterisation. K must be greater than 0."))
    end
    n = BioSequences.seq_data_len(A, K)
    if n !== N
        # This has been significantly changed conceptually from before. Now we
        # don't just check K, but *enforce* the most appropriate N for K.
        throw(ArgumentError("Bad kmer parameterisation. For K = $K, N should be $n"))
    end
end

@inline Base.length(x::Kmer{A,K,N}) where {A,K,N} = K
@inline Base.summary(x::Kmer{A,K,N}) where {A,K,N} = string(eltype(x), ' ', K, "-mer")

function Base.typemin(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}(unsafe, ntuple(i -> zero(UInt64), N))
end

function Base.typemax(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}((typemax(UInt64), ntuple(i -> typemax(UInt64), N - 1)...))
end

@inline function rand_kmer_data(::Type{Kmer{A,K,N}}, ::Val{true}) where {A,K,N}
    return Kmer{A,K,N}(ntuple(i -> rand(UInt64), Val{N}()))
end

@inline function rand_kmer_data(::Type{Kmer{A,K,N}}, ::Val{false}) where {A,K,N}
    ## All based on alphabet type of Kmer, so should constant fold.
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in 1:n_head
        bits = UInt64(BioSequences.encode(A(), rand(symbols(A()))))
        head = (head << bits_per_sym) | bits
    end
    # And the rest of the sequence
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            bits = UInt64(BioSequences.encode(A(), rand(symbols(A()))))
            body = (body << bits_per_sym) | bits
        end
        return body
    end
    return (head, tail...)
end

@inline function Base.rand(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    return Kmer{A,K,N}(rand_kmer_data(Kmer{A,K,N}, BioSequences.iscomplete(A())))
end

Base.rand(::Type{Kmer{A,K}}) where {A,K} = rand(kmertype(Kmer{A,K}))

function Base.rand(::Type{T}, size::Integer) where {T<:Kmer}
    return [rand(T) for _ in 1:size]
end

###
### Old Mer Base Functions - not transferred to new type.
###
#@inline encoded_data_type(::Type{Mer{A,K}}) where {A,K} = UInt64
#@inline encoded_data_type(::Type{BigMer{A,K}}) where {A,K} = UInt128
#@inline encoded_data_type(x::AbstractMer) = encoded_data_type(typeof(x))
#@inline encoded_data(x::AbstractMer) = reinterpret(encoded_data_type(typeof(x)), x)
#@inline ksize(::Type{T}) where {A,K,T<:AbstractMer{A,K}} = K
#@inline Base.unsigned(x::AbstractMer) = encoded_data(x)
#Base.:-(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) - y % encoded_data_type(x))
#Base.:+(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) + y % encoded_data_type(x))
#Base.:+(x::AbstractMer, y::AbstractMer) = y + x
#Alphabet(::Type{Mer{A,K} where A<:NucleicAcidAlphabet{2}}) where {K} = Any

include("indexing.jl")

#LongSequence{A}(x::Kmer{A,K,N}) where {A,K,N} = LongSequence{A}([nt for nt in x])
# Convenience method so as don't need to specify A in LongSequence{A}.
BioSequences.LongSequence(x::Kmer{A,K,N}) where {A,K,N} = LongSequence{A}(x)

include("predicates.jl")
include("counting.jl")
include("transformations.jl")

###
### Kmer de-bruijn neighbors
###

# TODO: Decide on this vs. old iterator pattern. I like the terseness of the code vs defining an iterator. Neither should allocate.
fw_neighbors(kmer::Kmer{A,K,N}) where {A<:DNAAlphabet,K,N} = ntuple(i -> pushlast(kmer, ACGT[i]), Val{4}())
fw_neighbors(kmer::Kmer{A,K,N}) where {A<:RNAAlphabet,K,N} = ntuple(i -> pushlast(kmer, ACGU[i]), Val{4}())
bw_neighbors(kmer::Kmer{A,K,N}) where {A<:DNAAlphabet,K,N} = ntuple(i -> pushfirst(kmer, ACGT[i]), Val{4}())
bw_neighbors(kmer::Kmer{A,K,N}) where {A<:RNAAlphabet,K,N} = ntuple(i -> pushfirst(kmer, ACGU[i]), Val{4}())

#=
# Neighbors on a de Bruijn graph
struct KmerNeighborIterator{S<:Kmer}
    x::S
end

"""
    neighbors(kmer::S) where {S<:Kmer}

Return an iterator through skip-mers neighboring `skipmer` on a de Bruijn graph.
"""
neighbors(kmer::Kmer) = KmerNeighborIterator{typeof(kmer)}(kmer)

Base.length(::KmerNeighborIterator) = 4
Base.eltype(::Type{KmerNeighborIterator{S}}) where {S<:Kmer} = S

function Base.iterate(it::KmerNeighborIterator{S}, i::UInt64 = 0) where {S<:Kmer}
    if i == 4
        return nothing
    else
        #return S((encoded_data(it.x) << 2) | i), i + 1
        return it.x << 1, i + one(UInt64)
    end
end
=#

###
### String literals
###

macro mer_str(seq, flag)
    seq′ = BioSequences.remove_newlines(seq)
    if flag == "dna" || flag == "d"
        T = kmertype(DNAKmer{length(seq′)})
        return T(seq′)
    elseif flag == "rna" || flag == "r"
        T = kmertype(RNAKmer{length(seq′)})
        return T(seq′)
    elseif flag == "aa" || flag == "a" || flag == "prot" || flag == "p" 
        T = kmertype(AAKmer{length(seq′)})
        return T(seq′)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

macro mer_str(seq)
    seq′ = BioSequences.remove_newlines(seq)
    T = kmertype(DNAKmer{length(seq′)})
    return T(seq′)
end

include("revtrans.jl")