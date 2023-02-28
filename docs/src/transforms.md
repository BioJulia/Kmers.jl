```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Indexing & modifying kmers

## Indexing

As `BioSequence` concrete subtypes, kmers can be indexed using integers

```jldoctest
julia> kmer = Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C)
DNA 5-mer:
TTAGC

julia> kmer[3]
DNA_A
```

You can also slice Kmers using UnitRanges:

```jldoctest
julia> kmer = Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C)
DNA 5-mer:
TTAGC

julia> kmer[1:3]
DNA 3-mer:
TTA
```

!!! warning
    Using slicing will introduce performance penalties in your code if
    you pass values of `i` that are not constants that can be propagated.

## Modifying sequences

Many modifying operations that are possible for some `BioSequences` such as
`LongSequence` are not possible for `Kmer`s, this is primarily due to the fact
`Kmer`s are an immutable struct.

However some non-mutating transformations are available:

```@docs
BioSequences.complement(::Kmer)
Base.reverse(::Kmer)
BioSequences.reverse_complement(::Kmer)
canonical
```