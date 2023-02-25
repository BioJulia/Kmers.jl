```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Iterating over kmers

When introducing the `Kmer` type we described kmers as contiguous sub-strings of
k nucleotides of some reference sequence.

This package therefore contains functionality for iterating over all the valid
`Kmers{A,K,N}` in a longer `BioSequence`.

```@docs
EveryKmer
EveryKmer{T,S}(seq::S, start::Int = firstindex(seq), stop::Int = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
EveryKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
SpacedKmers
SpacedKmers{T,S}(seq::S, step::Int, start::Int, stop::Int) where {T<:Kmer,S<:BioSequence}
SpacedKmers(seq::BioSequence{A}, ::Val{K}, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
EveryCanonicalKmer
EveryCanonicalKmer{T}(seq::S, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
EveryCanonicalKmer(seq::BioSequence{A}, ::Val{K}, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
SpacedCanonicalKmers
SpacedCanonicalKmers{T}(seq::S, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {T<:Kmer,S<:BioSequence}
SpacedCanonicalKmers(seq::BioSequence{A}, ::Val{K}, step::Int, start = firstindex(seq), stop = lastindex(seq)) where {A,K}
```


