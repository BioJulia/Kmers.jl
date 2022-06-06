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
SpacedKmers
EveryCanonicalKmer
SpacedCanonicalKmers
```


