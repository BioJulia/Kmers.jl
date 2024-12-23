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
Also, remember that `Kmer`s are only efficient when short (at most a few hundred symbols). Hence, entire exons or genes should probably not ever be represented by a `Kmer`, but rather as a `LongSequence` or `LongSubSeq` from BioSequences.jl.

### Reverse translation
Kmers.jl implements reverse translation, which maps an amino acid sequence to one or more RNA sequences.
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

Both functions take a genetic code as a keyword argument of the type `ReverseGeneticCode`. This object determines the mapping from amino acid to `CodonSet` - by default the [standard genetic code](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Standard_RNA_codon_table) is used - this mapping is used by nearly all organisms.

Only the reverse standard genetic code is defined in Kmers.jl.
To use another genetic code, build a `ReverseGeneticCode` object from an existing
`BioSequences.GeneticCode`:

```jldoctest
julia> code = BioSequences.pterobrachia_mitochondrial_genetic_code;

julia> rv_code = ReverseGeneticCode(code);

julia> seq = aa"KWLP";

julia> codonsets = reverse_translate(seq, rv_code)
4-element Vector{CodonSet}:
 CodonSet(0x0000000000000405)
 CodonSet(0x0500000000000000)
 CodonSet(0x50000000f0000000)
 CodonSet(0x0000000000f00000)

julia> codonsets == reverse_translate(seq) # default standard code
false
```

```@docs
ReverseGeneticCode
```
