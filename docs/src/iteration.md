```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```
## Iteration
Most applications of kmers extract multiple kmers from an underlying sequence.
To facilitate this, Kmers.jl implements a few basic kmer iterators, most of which are subtypes of `AbstractKmerIterator`.

The underlying sequence can be a `BioSequence`, `AbstractString`, or `AbstractVector{UInt8}`.
In the latter case, if the alphabet of the element type implements `BioSequences.AsciiAlphabet`, the vector will be treated as a vector of ASCII characters.

Similarly to the rules when constructing kmers directly, DNA and RNA is treated interchangeably when the underlying sequence is a `BioSequence`, but when the underlying sequence is a string or bytevector, `U` and `T` are considered different, and e.g. uracil cannot be constructed from a sequence containing `T`:

```jldoctest
julia> only(FwDNAMers{3}(rna"UGU"))
DNA 3-mer:
TGT

julia> only(FwDNAMers{3}("UGU"))
ERROR:
[...]
```

The following kmer iterators are implemented:

### `FwKmers`
The most basic kmer iterator is `FwKmers`, which simply iterates every kmer, in order:

```@docs
FwKmers
FwDNAMers
FwRNAMers
FwAAMers
```

### `FwRvIterator`
This iterates over a nucleic acid sequence. For every kmer it encounters, it outputs the kmer and its reverse complement.

```@docs
FwRvIterator
```

### `CanonicalKmers`
This iterator is similar to [`FwKmers`](@ref), however, for each `Kmer` encountered, it returns the _canonical_ kmer.

The canonical kmer is defined as the lexographically smaller of a kmer and its reverse complement.
That is, if [`FwKmers`](@ref) would iterate `TCAC`, then [`CanonicalKmers`](@ref) would return `GTGA`, as this is the reverse complement of `TCAC`, and is before `TCAC` in the alphabet.

[`CanonicalKmers`](@ref) is useful for summarizing the kmer composition of sequences whose strandedness is unknown.

```@docs
CanonicalKmers
CanonicalDNAMers
CanonicalRNAMers
```

### `UnambiguousKmers`
[`UnambiguousKmers`](@ref) iterates unambiguous nucleotides (that is, kmers of the alphabets `DNAAlphabet{2}` or `RNAAlphabet{2}`).
Any kmers containing [ambiguous nucleotides](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC341218/) such as `W` or `N` are skipped.

```@docs
UnambiguousKmers
UnambiguousDNAMers
UnambiguousRNAMers
```

### `SpacedKmers`
The [`SpacedKmers`](@ref) iterator iterates kmers with a fixed step size between k-mers.
For example, for a K of 4, and a step size of 3, the output kmers would overlap with a single nucleotide, like so:

```
seq: TGATGCGTAGTG
     TGCT
        TGCG
           GTAG
```

Hence, if `FwKmers` are analogous to `UnitRange`, `SpacedKmers` is analogous to `StepRange`.

```@docs
SpacedKmers
SpacedDNAMers
SpacedRNAMers
SpacedAAMers
```

The convenience functions [`each_codon`](@ref) return `SpacedKmers` with a K value of 3 and step size of 3:

```@docs
each_codon
```

## The `AbstractKmerIterator` interface
It's very likely that users of Kmers.jl need to implement their own custom kmer iterators, in which case they should subtype [`AbstractKmerIterator`](@ref).

```@docs
AbstractKmerIterator
```

At the moment, there is no real interface implemented for this abstract type,
other than that `AbstractKmerIterator{A, K}` needs to iterate `Kmer{A, K}`.
