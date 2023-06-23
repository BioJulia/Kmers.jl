```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```
## FAQ
### Why can kmers not be compared to biosequences?
It may be surprising that kmers cannot be compared to other biosequences:

```jldoctest
julia> dna"TAG" == mer"TAG"d
ERROR: MethodError
[...]
```

In fact, this is implemented by a manually thrown `MethodError`; the generic case `Base.:==(::BioSequence, ::BioSequence)` is defined.

The reason for this is the consequence of the following limitations:
* `isequal(x, y)` implies `hash(x) == hash(y)`
* `isequal(x, y)` and `x == y` ought to be identical for well-defined elements (i.e. in the absence of `missing`s and `NaN`s etc.)
* `hash(::Kmer)` must be absolutely maximally efficient

If kmers were to be comparable to `BioSequence`, then the hashing of `BioSequence` should follow `Kmer`, which practically speaking would mean that all biosequences would need to be recoded to `Kmer`s before hashing.

### Why isn't there an iterator of unambiguous, canonical kmers or spaced, canonical kmers?
Any iterator of nucleotide kmers can be made into a canonical kmer iterator by simply calling `canonical` on its output kers.

The `CanonicalKmers` iterator is special cased, because with a step size of 1, it is generally faster to build the next kmer by storing both the reverse and forward kmer, then creating the next kmer by prepending/append the next symbol.

However, with a larger step size, it becomes more efficient to build the forward kmer, then reverse-complement the whole kmer.

### Why isn't there an iterator of skipmers/minimizers/k-min-mers, etc?
The concept of kmers have turned out to be remarkably flexible and useful in bioinformatics, and have spawned a neverending stream of variations.
We simply can't implement them all.

However, see the section [Building kmer replacements](@ref replacements) on how to implement them
as a user of Kmers.jl yourself.
