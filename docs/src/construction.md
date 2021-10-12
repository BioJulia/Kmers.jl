```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Construction & conversion

Here we will showcase the various ways you can construct Kmers.
Because `Kmer` has type parameters for their length, some methods of construction
are type stable, and others - notable ones that get the length at runtime, are not.

Each section below is split into type stable and non type stable examples so you
can learn the difference.

## From strings

Kmers can be constructed from strings using their constructors:

### Type stable

```jldoctest
julia> DNAKmer{8}("ATCGATCG")
DNA 8-mer:
ATCGATCG

julia> RNAKmer{8}("AUCGAUCG")
RNA 8-mer:
AUCGAUCG
```

### Not type stable

```jldoctest
julia> DNAKmer("ATCGATCG")
DNA 8-mer:
ATCGATCG

julia> RNAKmer("AUCGAUCG")
RNA 8-mer:
AUCGAUCG
```

!!! tip
    The difference here between achieving type stability or not is in
    providing the length (if you know it) as a parameter

## From arrays of BioSymbols

Kmers can be constructed using vectors or arrays of a BioSymbol type:

### Type stable

```jldoctest
julia> DNAKmer{5}([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
DNA 5-mer:
TTAGC

julia> RNAKmer{5}([RNA_U, RNA_U, RNA_A, RNA_G, RNA_C])
RNA 5-mer:
UUAGC
```

### Not type stable

```jldoctest
julia> DNAKmer([DNA_T, DNA_T, DNA_A, DNA_G, DNA_C])
DNA 5-mer:
TTAGC

julia> RNAKmer([RNA_U, RNA_U, RNA_A, RNA_G, RNA_C])
RNA 5-mer:
UUAGC
```

!!! tip
    The difference here between achieving type stability or not is in
    providing the length (if you know it) as a parameter


## From other sequences

### Type stable

```jldoctest
julia> seq = LongSequence("TTAGC")
5nt DNA Sequence:
TTAGC

julia> DNAKmer{5}(seq)
DNA 5-mer:
TTAGC
```

### Not type stable

```jldoctest
julia> seq = LongSequence("TTAGC")
5nt DNA Sequence:
TTAGC

julia> DNAKmer(seq)
DNA 5-mer:
TTAGC
```

!!! tip
    The difference here between achieving type stability or not is in
    providing the length (if you know it) as a parameter

