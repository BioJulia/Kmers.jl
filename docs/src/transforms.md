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
julia> seq = Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C)
DNA 5-mer:
TTAGC

julia> seq[3]
DNA_A
```

Currently, indexing Kmers using arbitrary ranges is not implemented because it
is not possible to do in a type-stable way.

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