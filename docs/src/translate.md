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

## Reverse translation
Reverse translation (or "revtrans", for short) refers to the mapping from amino acids back to the set of RNA codons that code for the given amino acid, under a given genetic code.
There is no known natural process of revtrans, but it can be useful to do _in silico_.

In Kmers.jl, revtrans is done through the `reverse_translate` function.
This takes an amino acid sequence and produces a `Vector{CodonSet}`, where `CodonSet <: AbstractSet{RNACodon}`.
Alternatively, it takes an amino acid and produces a `CodonSet`.

A reverse genetic code can optionally be specified as the second argument.
If not provided, it default to the reverse standard genetic code.

### Example of reverse translation
```julia
julia> reverse_translate(AA_W) # default to standard genetic code
Kmers.CodonSet with 1 element:
  UGG

julia> code = ReverseGeneticCode(BioSequences.trematode_mitochondrial_genetic_code);

julia> reverse_translate(AA_W, code)
Kmers.CodonSet with 2 elements:
  UGA
  UGG
```

### Important notes on reverse translation
* `AA_Gap` cannot be reverse translated. Attempting so throws an error
* In cells, `AA_O` and `AA_U` are encoded by dynamic overloading of the codons `UAG` and `UGA`, respectively.
  Because these codons normally code for `AA_Term`, the forward genetic code returns `AA_Term` for these codons.
  However, we can unambiguously reverse translate them, so these amino acids translate to codonsets with these
  precise codons.
* Ambiguous amino acids translate to the union of the possible amino acids. For example, if `AA_L` translate to set `S1`,
  and `AA_I` translate to `S2`, then `AA_J` translate to `union(S1, S2)`.

```@docs
Kmers.CodonSet
Kmers.ReverseGeneticCoTranslatingde
reverse_translate
reverse_translate!
```
