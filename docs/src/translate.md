```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Translating and reverse translating

## Translating
Just like other `BioSequence`s, `Kmer`s of RNA or DNA alphabets can be efficiently translated to amino acids:

```
julia> kmer = RNAKmer("AUGGGCCACUGA");

julia> translate(kmer)
AminoAcid 4-mer:
MGH*
```

For more information on translation and different genetic codes, see the documentation of BioSequences.jl.
