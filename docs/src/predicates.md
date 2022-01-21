```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Predicates

The following predicate functions from BioSequences.jl are compatible with `Kmer`s.
Some have an optimised method defined in Kmers.jl.

```@docs
isrepetitive
ispalindromic
hasambiguity
iscanonical
```