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

###
### _build_kmer_data
###

#=
These are (hopefully!) very optimised kernel functions for building kmer internal
data from individual elements or from sequences. Kmers themselves are static,
tuple-based structs, and so I really didn't want these functions to create memory
allocations or GC activity through use of vectors an such, for what should be
the creation of a single, rather simple value.
=#

"""
    _build_kmer_data(::Type{Kmer{A,K,N}}, seq::LongSequence{A}, from::Int = 1) where {A,K,N}

Construct a ntuple of the bits data for an instance of a Kmer{A,K,N}.

This particular method is specialised for LongSequences, and for when the Kmer
and LongSequence types used, share the same alphabet, since a lot of encoding /
decoding can be skipped, and the problem is mostly one of shunting bits around.
"""
@inline function _build_kmer_data(::Type{Kmer{A,K,N}}, seq::LongSequence{A}, from::Int = 1) where {A,K,N}
    checkmer(Kmer{A,K,N})
    
    bits_per_sym = BioSequences.bits_per_symbol(A()) # Based on alphabet type, should constant fold.
    n_head = elements_in_head(Kmer{A,K,N}) # Based on kmer type, should constant fold.
    n_per_chunk = per_word_capacity(Kmer{A,K,N}) # Based on kmer type, should constant fold.
    
    if from + K - 1 > length(seq)
        return nothing
    end
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in from:(from + n_head - 1)
        bits = UInt64(BioSequences.extract_encoded_element(seq, i))
        head = (head << bits_per_sym) | bits
    end
    
    # And the rest of the sequence
    idx = Ref(from + n_head)
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            bits = UInt64(BioSequences.extract_encoded_element(seq, idx[]))
            body = (body << bits_per_sym) | bits
            idx[] += 1
        end
        return body
    end
    
    # Put head and tail together
    return (head, tail...)
end



###
### Constructors
###

"""
    Kmer{A,K,N}(itr) where {A,K,N}

Construct a `Kmer{A,K,N}` from an iterable.

The most generic constructor.

Currently the iterable must have `length` & support `getindex` with integers.

# Examples

```jldoctest
julia> ntseq = LongSequence("TTAGC") # 4-bit DNA alphabet
5nt DNA Sequence:
TTAGC

julia> DNAKmer{5}(ntseq) # 2-Bit DNA alphabet
DNA 5-mer:
TTAGC
```
"""
function Kmer{A,K,N}(itr) where {A,K,N}
    checkmer(Kmer{A,K,N})
    
    seqlen = length(itr)
    if seqlen != K
        throw(ArgumentError("itr does not contain enough elements ($seqlen ≠ $K)"))
    end
    
    ## All based on alphabet type of Kmer, so should constant fold.
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in 1:n_head
        sym = convert(eltype(Kmer{A,K,N}), itr[i])
        # Encode will throw if it cant encode an element.
        head = (head << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
    end
    
    # And the rest of the sequence
    idx = Ref(n_head + 1)
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for i in 1:n_per_chunk
            sym = convert(eltype(Kmer{A,K,N}), itr[idx[]])
            # Encode will throw  if it cant encode an element.
            body = (body << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
            idx[] += 1
        end
        return body
    end
    
    data = (head, tail...)
    
    return Kmer{A,K,N}(data)
end

"""
    Kmer{A,K,N}(seq::BioSequence{A})

Construct a `Kmer{A,K,N}` from a `BioSequence{A}`.

This particular method is specialised for BioSequences, and for when the Kmer
and BioSequence types used, share the same alphabet, since a lot of encoding /
decoding can be skipped, and the problem is mostly one of shunting bits around.
In the case where the alphabet of the Kmer and the alphabet of the BioSequence
differ, dispatch to the more generic constructor occurs instead.

# Examples

```jldoctest
julia> ntseq = LongSequence{DNAAlphabet{2}}("TTAGC") # 2-bit DNA alphabet
5nt DNA Sequence:
TTAGC

julia> DNAKmer{5}(ntseq) # 2-Bit DNA alphabet
DNA 5-mer:
TTAGC
```
"""
@inline function Kmer{A,K,N}(seq::BioSequence{A}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    
    seqlen = length(seq)
    if seqlen != K
        throw(ArgumentError("seq is not the correct length ($seqlen ≠ $K)"))
    end
    
    ## All based on alphabet type of Kmer, so should constant fold.
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in 1:n_head
        bits = UInt64(BioSequences.extract_encoded_element(seq, i))
        head = (head << bits_per_sym) | bits
    end
    
    # And the rest of the sequence
    idx = Ref(n_head + 1)
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            bits = UInt64(BioSequences.extract_encoded_element(seq, idx[]))
            body = (body << bits_per_sym) | bits
            idx[] += 1
        end
        return body
    end
    
    data = (head, tail...)
    
    return Kmer{A,K,N}(data)
end


# Convenience version of function above so you don't have to work out correct N.
"""
    Kmer{A,K}(itr) where {A,K}

Construct a `Kmer{A,K,N}` from an iterable.

This is a convenience method which will work out the correct `N` parameter, for
your given choice of `A` & `K`.
"""
@inline function Kmer{A,K}(itr) where {A,K}
    T = kmertype(Kmer{A,K})
    return T(itr)
end

"""
    Kmer{A}(itr) where {A}

Construct a `Kmer{A,K,N}` from an iterable.

This is a convenience method which will work out K from the length of `itr`, and
the correct `N` parameter, for your given choice of `A` & `K`.

!!! warning
    Since this gets K from runtime values, this is gonna be slow!
"""
@inline Kmer{A}(itr) where {A} = Kmer{A,length(itr)}(itr)
@inline Kmer(seq::BioSequence{A}) where A = Kmer{A}(seq)

function Kmer{A1}(seq::BioSequence{A2}) where {A1 <: NucleicAcidAlphabet, A2 <: NucleicAcidAlphabet}
    kmertype(Kmer{A1, length(seq)})(seq)
end

@inline function Kmer{A}(nts::Vararg{Union{DNA, RNA}, K}) where {A <: NucleicAcidAlphabet, K}
    return kmertype(Kmer{A, K})(nts)
end

"""
    Kmer(nts::Vararg{DNA,K}) where {K}

Construct a Kmer from a variable number `K` of DNA nucleotides.

# Examples

```jldoctest
julia> Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C)
DNA 5-mer:
TTAGC
```
"""
@inline Kmer(nt::DNA, nts::Vararg{DNA}) = DNAKmer((nt, nts...))

"""
    Kmer(nts::Vararg{RNA,K}) where {K}

Construct a Kmer from a variable number `K` of RNA nucleotides.

# Examples

```jldoctest
julia> Kmer(RNA_U, RNA_U, RNA_A, RNA_G, RNA_C)
DNA 5-mer:
UUAGC
```
"""
@inline Kmer(nt::RNA, nts::Vararg{RNA}) = RNAKmer((nt, nts...))


"""
    Kmer(seq::String)

Construct a DNA or RNA kmer from a string.

!!! warning
    As a convenience method, this derives the `K`, `Alphabet`, and `N` parameters
    for the `Kmer{A,K,N}` type from the input string.

# Examples

```jldoctest
julia> Kmer("TTAGC")
DNA 5-mer:
TTAGC
```
"""
@inline function Kmer(seq::String)
    seq′ = BioSequences.remove_newlines(seq)
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
    return Kmer{A,K,BioSequences.seq_data_len(A, K)}
end
@inline kmertype(::Type{Kmer{A,K,N}}) where {A,K,N} = Kmer{A,K,N}

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


"""
    Base.rand(::Type{Kmer{A,K,N}}) where {A,K,N}
    Base.rand(::Type{Kmer{A,K}}) where {A,K}

Create a random kmer of a specified alphabet and length

# Examples
```julia
julia> rand(Kmer{DNAAlphabet{2}, 3})
BioSymbols.DNA 3-mer:
ACT

```
"""
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