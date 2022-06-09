```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Kmer types

Bioinformatic analyses make extensive use of kmers.
Kmers are contiguous sub-strings of k nucleotides of some reference sequence. 

They are used extensively in bioinformatic analyses as an informational unit.
This concept popularised by short read assemblers. 
Analyses within the kmer space benefit from a simple formulation of the sampling
problem and direct in-hash comparisons.

```@docs
Kmer
```

The following aliases are also defined:

```@docs
DNAKmer
DNA27mer
DNA31mer
DNA63mer
RNAKmer
RNA27mer
RNA31mer
RNA63mer
```
