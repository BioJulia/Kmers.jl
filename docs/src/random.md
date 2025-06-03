```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
    using Random
end
```
# Random genration of Kmers
Kmers.jl supports efficient random generation of kmers through a package extension, when the `Random` stdlib is loaded.

## Sampling a symbol from a kmer
```jldoctest; filter = r"^AA_[PLQWHY]"
julia> m = mer"PLQWHY"a;

julia> rand(m)
AA_L

julia> rand(m)
AA_H
```

# Sampling a random kmer type
When sampling a random kmer from a type, kmers are sampled from the following distributions:

* For `T <: AAKmer`, the kmer symbols is uniformly sampled from the 20 standard proteogenic amino acids `ACDEFGHIKLMNPQRSTVWY`. This excludes amino acids such as selenocysteine, and pseudo-aas such as `AA_Term` and `AA_Gap`
* For `T <: DNA{<:NucleicAcidAlphabet{4}}`, the kmer symbols is uniformly sampled from the 4 unambiguous nucleotides. E.g. for `T <: Kmer{DNAAlphabet{4}}`, the only symbols are `A`, `C`, `G` and `T`.
* For all other kmer types, the symbols are sampled uniformly from all symbols in its alphabet.

You can sample from either a concrete type, or a type with the trailing `N` type parameter elided:
```jldoctest; filter = r"[A-Z]+"s
julia> rand(DNAKmer{9})
DNA 9-mer:
TGGGAGCCA

julia> rand(AAKmer{4, 1})
AminoAcid 4-mer:
RRIP
```