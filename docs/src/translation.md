```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```

## Translation
`Kmer`s can be translated using the `translate` function exported by `BioSequences`:

```jldoctest
julia> translate(mer"UGCUUGAUC"r)
AminoAcid 3-mer:
CLI
```

Since `Kmer`s are immutable, the in-place `translate!` function is not implemented for `Kmers`.
Also, remember that `Kmer`s are only efficient when short (at most a few hundred symbols). Hence, entire exons or genes should probably be represented by `LongSequence` or `LongSubSeq`.

### Reverse translation
Kmers.jl implements reverse translation, in which an amino acid sequence is translated to an RNA sequence.
While this process doesn't occur naturally (as far as we know), it is still useful for some analyses.

Since genetic codes are degenerate, i.e. multiple codons code for the same amino acid, reverse translating a sequence does not return a nucleic acid sequence, but a vector of `CodonSet`:

```@docs
reverse_translate
CodonSet
```

`CodonSet` is an efficiently implemented `AbstractSet{RNACodon}` (and remember, `RNACodon` is an alias for `RNAKmer{3, 1}`).

To avoid allocating a new `Vector`, you can use `reverse_translate!`:

```@docs
reverse_translate!
```

Both functions take a genetic code as a keyword argument of the type `ReverseGeneticCode`. This object determines the mapping from amino acid to `CodonSet` - by default the [standard genetic code](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Standard_RNA_codon_table) is used - this mapping is used by nearly all organisms:

```@docs
ReverseGeneticCode
```