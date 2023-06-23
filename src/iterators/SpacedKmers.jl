"""
    SpacedKmers{A <: Alphabet, K, J, S} <: AbstractKmerIterator{A, K}

Iterator of kmers with step size. `J` signifies the step size, `S`
the type of the underlying sequence, and the eltype of the iterator
is `Kmer{A, K, N}` with the appropriate `N`.

For example, a `SpacedKmers{AminoAcidAlphabet, 3, 5, Vector{UInt8}}` sampling
over `seq::Vector{UInt8}` will sample all kmers corresponding to
`seq[1:3], seq[6:8], seq[11:13]` etc.

See also: [`each_codon`](@ref), [`FwKmers`](@ref)

# Examples:
```jldoctest
julia> collect(SpacedDNAMers{3, 2}("AGCGTATA"))
3-element Vector{Kmer{DNAAlphabet{2}, 3, 1}}:
 AGC
 CGT
 TAT
```
"""
struct SpacedKmers{A <: Alphabet, K, J, S} <: AbstractKmerIterator{A, K}
    seq::S

    function SpacedKmers{A, K, J, S}(seq::S) where {A, K, J, S}
        K isa Int || error("K must be an Int")
        K > 0 || error("K must be at least 1")
        J isa Int || error("J must be an Int")
        J > 0 || error("J must be at least 1")
        new{A, K, J, S}(seq)
    end
end

source_type(::Type{SpacedKmers{A, K, J, S}}) where {A, K, J, S} = S
stepsize(::SpacedKmers{A, K, J}) where {A, K, J} = J

@inline function Base.length(it::SpacedKmers{A, K, J}) where {A, K, J}
    src = used_source(RecodingScheme(A(), source_type(typeof(it))), it.seq)
    L = length(src)
    L < K ? 0 : div((L - K), J) + 1
end

SpacedKmers{A, K, J}(s) where {A <: Alphabet, K, J} = SpacedKmers{A, K, J, typeof(s)}(s)

"`SpacedDNAMers{K, J, S}`: Alias for `SpacedKmers{DNAAlphabet{2}, K, J, S}`"
const SpacedDNAMers{K, J, S} = SpacedKmers{DNAAlphabet{2}, K, J, S}

"`SpacedRNAMers{K, J, S}`: Alias for `SpacedKmers{RNAAlphabet{2}, K, J, S}`"
const SpacedRNAMers{K, J, S} = SpacedKmers{RNAAlphabet{2}, K, J, S}

"`SpacedAAMers{K, J, S}`: Alias for `SpacedKmers{AminoAcidAlphabet, K, J, S}`"
const SpacedAAMers{K, J, S} = SpacedKmers{AminoAcidAlphabet, K, J, S}

"""
    each_codon(s::BioSequence{<:Union{DNAAlphabet, RNAAlphabet}})
    each_codon(::Type{<:Union{DNA, RNA}}, s)

Construct an iterator of nucleotide 3-mers with step size 3 from `s`.
The sequence `s` may be an RNA or DNA biosequence, in which case the element
type is inferred, or the element type may be specified explicitly, in which
case `s` may be a byte-like sequence such as a `String` or `Vector{UInt8}`.

This function returns [`SpacedKmers`](@ref) iterator.

See also: [`SpacedKmers`](@ref)

Examples:
```jldoctest
julia> collect(each_codon(DNA, "TGACGATCGAC"))
3-element Vector{Kmer{DNAAlphabet{2}, 3, 1}}:
 TGA
 CGA
 TCG
```
"""
each_codon(::Type{DNA}, s) = SpacedDNAMers{3, 3}(s)
each_codon(::Type{RNA}, s) = SpacedRNAMers{3, 3}(s)

each_codon(s::BioSequence{<:DNAAlphabet}) = SpacedDNAMers{3, 3}(s)
each_codon(s::BioSequence{<:RNAAlphabet}) = SpacedRNAMers{3, 3}(s)

@inline function Base.iterate(it::SpacedKmers{A}, state...) where {A}
    iterate_kmer(RecodingScheme(A(), source_type(typeof(it))), it, state...)
end

# TODO: Maybe in all kmer iterators, instantiate it with the source type,
# so we don't have to get the source type in functions (and thus
# it is allwoed to be a costly operation).
# However, this means we instantiate e.g. a FwKmers{A, K, S} and change S
# in the source type in the constructor
@inline function iterate_kmer(
    R::RecodingScheme,
    it::SpacedKmers{A, K},
) where {A <: Alphabet, K}
    length(it.seq) < K && return nothing
    kmer = unsafe_extract(
        R,
        eltype(it),
        used_source(RecodingScheme(A(), source_type(typeof(it))), it.seq),
        1,
    )
    next_index = 1 + max(stepsize(it), K)
    (kmer, (kmer, next_index))
end

# Here, we need to convert to an abstractvector
# TODO: This function and the one above can be merged with the FwKmers one?
@inline function iterate_kmer(
    R::AsciiEncode,
    it::SpacedKmers{A, K, J, S},
) where {A <: Alphabet, K, J, S}
    src = used_source(RecodingScheme(A(), S), it.seq)
    Base.require_one_based_indexing(src)
    length(src) < K && return nothing
    kmer = unsafe_extract(R, eltype(it), src, 1)
    next_index = 1 + max(stepsize(it), K)
    (kmer, (kmer, next_index))
end

@inline function iterate_kmer(
    ::RecodingScheme,
    it::SpacedKmers{A, K, J, S},
    state,
) where {A, K, S, J}
    src = used_source(RecodingScheme(A(), S), it.seq)
    R = RecodingScheme(A(), S)
    Base.require_one_based_indexing(src)
    (kmer, i) = state
    i > lastindex(src) - min(K, J) + 1 && return nothing
    next_i = i + J
    # This branch should be resolved statically
    if J ≥ K
        kmer = unsafe_extract(R, eltype(it), src, i)
    else
        kmer = unsafe_shift_from(R, kmer, src, i, Val{J}())
    end
    (kmer, (kmer, next_i))
end
