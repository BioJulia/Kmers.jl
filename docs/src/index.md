## Kmers.jl
Kmers.jl provide the `Kmer <: BioSequence` type which implement the concept of a
[k-mer](https://en.wikipedia.org/wiki/K-mer), a biological sequence of exactly length `k`.


K-mers are used frequently in bioinformatics because, when k is small and known at
compile time, these sequences can be efficiently represented as integers and stored
directly in CPU registers, allowing for much more efficient computation than arbitrary-length sequences.

In Kmers.jl, the `Kmer` type is psrameterized by its length, and its data is stored in an `NTuple`. This makes `Kmers` bitstypes and highly efficient.

Conceptually, one may use the following analogy:
* `BioSequence` is like `AbstractVector`
* `LongSequence` is like `Vector`
* `Kmer` is like [`SVector`](https://github.com/JuliaArrays/StaticArrays.jl) from `StaticArrays`

Kmers.jl is tightly coupled to the
[`BioSequences.jl`](https://github.com/BioJulia/BioSequences.jl) package,
and relies on its internals.
Hence, you should expect strict compat bounds on BioSequences.jl.

## Usage
### ⚠️ WARNING ⚠️
`Kmer`s are parameterized by their length. That means any operation on `Kmer`s that change their length, such as `push`, `pop`, slicing, or masking (logical indexing) will be **type unstable** and hence slow and memory inefficient, unless you write your code in such as way that the compiler can use constant folding.

Kmers.jl is intended for high-performance computing. If you do not need the extra performance that register-stored sequences provide, you might consider using `LongSequence` from BioSequences.jl instead

## Installation
You can install BioSequences from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```julia
pkg> add Kmers
```

If you are interested in the cutting edge of development, please check out
the master branch to try new features before release.

## Contributing
We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.

## Questions?
If you have a question about contributing or using BioJulia software, come
on over and chat to us on [the Julia Slack workspace](https://julialang.org/slack/), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).
