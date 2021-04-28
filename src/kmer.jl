###
### Kmer Type definition
###

# Include some basic tuple bitflipping ops - the secret sauce to efficiently
# manipping Kmer's static data. 
include("tuple_bitflipping.jl")

"""
    Kmer{A<:Alphabet,K,N} <: BioSequence{A}

A parametric, immutable, bitstype for representing Kmers - short sequences.
Given the number of Kmers generated from raw sequencing reads, avoiding
repetetive memory allocation and triggering of garbage collection is important,
as is the ability to effectively pack Kmers into arrays and similar collections.
In julia this means an immutable bitstype must represent such shorter Kmer
sequences. Thankfully this is not much of a limitation - kmers are rarely
manipulated and so by and large don't have to be mutable like `LongSequence`s.
Excepting their immutability, they fulfill the rest of the API and behaviours
expected from a concrete `BioSequence` type, and non-mutating transformations
of the type are still defined.

!!! warning
    Given their immutability, `setindex` and mutating sequence transformations
    are not implemented for kmers e.g. `reverse_complement!`. 
!!! tip
    Note that some sequence transformations that are not mutating are
    available, since they can return a new kmer value as a result e.g.
    `reverse_complement`. 
"""
struct Kmer{A<:Alphabet,K,N} <: BioSequence{A}
    data::NTuple{N,UInt64}
    
    function Kmer{A,K,N}(data::NTuple{N,UInt64}) where {A<:Alphabet,K,N}
        checkmer(Kmer{A,K,N})
        # TODO: Decide on whether this method should always mask the (64N - 2K)
        # MSBs of the input tuple, as we do that in quite a few cases before
        # calling this constructor: see the typemin, typemax, rand, and transformations.jl
        return new(_cliphead(n_unused(Kmer{A,K,N}) * bits_per_symbol(A()), data...))
    end
end

# Create a blank ntuple of appropriate length for a given Kmer with N.
@inline blank_ntuple(::Type{Kmer{A,K,N}}) where {A,K,N} = ntuple(x -> zero(UInt64), Val{N}())


###
### Constructors
###

# Create a Mer from a sequence.
function Kmer{A,K,N}(nucs) where {A,K,N}
    seqlen = length(nucs)
    if seqlen != K
        throw(ArgumentError("nucs does not contain the correct number of nucleotides ($seqlen ≠ $K)"))
    end
    data = _build_kmer_data(nucs, Kmer{A,K,N})
    return Kmer{A,K,N}(data)
end

# Convenience version of function above so you don't have to work out correct N.
function Kmer{A,K}(nucs) where {A,K}
    T = kmertype(Kmer{A,K})
    return T(nucs)
end

# Convenience version of function above so you don't have to provide K.
# BUT THIS WILL GET K FROM RUNTIME VALUES - YA CODE IS GONNA BE SLLLOOOOOOWWWW.
Kmer{A}(nucs) where {A,K} = Kmer{A,length(nucs)}(nucs)

@inline function _build_kmer_data(nucs, ::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    # Construct the head.
    elements_per_chunk = div(64, bits_per_symbol(A()))
    bases_in_head = div(64 - (64N - (bits_per_symbol(A()) * K)), bits_per_symbol(A()))
    head = zero(UInt64)
    @inbounds for i in 1:bases_in_head
        nt = convert(eltype(Kmer{A,K,N}), nucs[i])
        if isgap(nt)
            throw(ArgumentError("cannot create a mer with gaps"))
        end
        head = (head << bits_per_symbol(A())) | UInt64(encode(A(), nt))
    end
    # And the rest of the sequence
    idx = Ref(bases_in_head + 1)
    
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for i in 1:elements_per_chunk
            nt = convert(eltype(Kmer{A,K,N}), nucs[idx[]])
            if isgap(nt)
                throw(ArgumentError("cannot create a mer with gaps"))
            end
            body = (body << bits_per_symbol(A())) | UInt64(encode(A(), nt))
            idx[] += 1
        end
        return body
    end
    
    return (head, tail...)
end

@inline function _build_kmer_data2(::Type{Kmer{A,K,N}}, nucs) where {A,K,N}
    checkmer(Kmer{A,K,N})
    kmer_data = blank_ntuple(Kmer{A,K,N})
    stop = length(nucs)
    while i ≤ stop
        nt = convert(eltype(A), @inbounds nucs[i])
        bits = encode(A(), nt)
        kmer_data = leftshift_carry(kmer_data, 2, fbits)
        filled += 1
        if filled == K
            return kmer_data
        end
        i += 1
    end
    return nothing
end











@inline function Kmer(nts::Vararg{DNA,K}) where {K}
    T = kmertype(DNAKmer{K})
    data = _build_kmer_data(nts, T)
    return T(data)
end

@inline function Kmer(nts::Vararg{RNA,K}) where {K}
    T = kmertype(RNAKmer{K})
    data = _build_kmer_data(nts, T)
    return T(data)
end

@inline function Kmer(seq::String)
    seq′ = remove_newlines(seq)
    hast = false
    hasu = false
    for c in seq′
        hast |= ((c == 'T') | (c == 't'))
        hasu |= ((c == 'U') | (c == 'u'))
    end
    if (hast & hasu) | (!hast & !hasu)
        throw(ArgumentError("Can't detect alphabet type from string"))
    end
    A = ifelse(hast & !hasu, DNAAlphabet{2}, RNAAlphabet{2})
    return Kmer{A,length(seq′)}(seq′)
end


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
    return Kmer{A,K,seq_data_len(A, K)}
end
@inline kmertype(::Type{Kmer{A,K,N}}) where {A,K,N} = Kmer{A,K,N}

# Aliases
"Shortcut for the type `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K,N} = Kmer{DNAAlphabet{2},K,N}

"Shortcut for the type `DNAKmer{27,1}`"
const DNAKmer27 = DNAKmer{27,1}

"Shortcut for the type `DNAKmer{31,1}`"
const DNAKmer31 = DNAKmer{31,1}

"Shortcut for the type `DNAKmer{63,2}`"
const DNAKmer63 = DNAKmer{63,2}

"Shortcut for the type `Kmer{RNAAlphabet{2},K,N}`"
const RNAKmer{K,N} = Kmer{RNAAlphabet{2},K,N}

"Shortcut for the type `RNAKmer{27,1}`"
const RNAKmer27 = RNAKmer{27,1}

"Shortcut for the type `RNAKmer{31,1}`"
const RNAKmer31 = RNAKmer{31,1}

"Shortcut for the type `RNAKmer{63,2}`"
const RNAKmer63 = RNAKmer{63,2}

"Shortcut for the type `Kmer{AminoAcidAlphabet,K,N}`"
const AAKmer{K,N} = Kmer{AminoAcidAlphabet,K,N}

const DNACodon = DNAKmer{3,1}
const RNACodon = RNAKmer{3,1}

###
### Base Functions
###

@inline ksize(::Type{Kmer{A,K,N}}) where {A,K,N} = K
@inline capacity(::Type{Kmer{A,K,N}}) where {A,K,N} = div(64N, bits_per_symbol(A()))
@inline capacity(seq::Kmer) = capacity(typeof(seq))
@inline n_unused(::Type{Kmer{A,K,N}}) where {A,K,N} = capacity(Kmer{A,K,N}) - K
@inline n_unused(seq::Kmer) = n_unused(typeof(seq))

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
    n = seq_data_len(A, K)
    if n !== N
        # This has been significantly changed conceptually from before. Now we
        # don't just check K, but *enforce* the most appropriate N for K.
        throw(ArgumentError("Bad kmer parameterisation. For K = $K, N should be $n"))
    end
end

@inline Base.length(x::Kmer{A,K,N}) where {A,K,N} = K
@inline Base.summary(x::Kmer{A,K,N}) where {A,K,N} = string(eltype(x), ' ', K, "-mer")

function Base.typemin(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}(ntuple(i -> zero(UInt64), N))
end

function Base.typemax(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}((typemax(UInt64), ntuple(i -> typemax(UInt64), N - 1)...))
end

function Base.rand(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}((rand(UInt64), ntuple(i -> rand(UInt64), N - 1)...))
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
LongSequence(x::Kmer{A,K,N}) where {A,K,N} = LongSequence{A}(x)

include("predicates.jl")
include("counting.jl")
include("transformations.jl")


###
### Kmer de-bruijn neighbors
###

# TODO: Decide on this vs. old iterator pattern. I like the terseness of the code vs defining an iterator. Neither allocate.
fw_neighbors(kmer::Kmer) = ntuple(i -> pushlast(kmer, ACGT[i]), Val{4}())
bw_neighbors(kmer::Kmer) = ntuple(i -> pushfirst(kmer, ACGT[i]), Val{4}())

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
    seq′ = remove_newlines(seq)
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
    seq′ = remove_newlines(seq)
    T = kmertype(DNAKmer{length(seq′)})
    return T(seq′)
end
