```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```

## Dynamically sized kmers
Sometimes, one requires using kmers of differing lengths in the same workload.
An example could be representing primers, which can be of length 18-24.
Here, using the `Kmer` type would cause code to specialize on each kmer length.
Beside causing both excessive compilation and code generation, it will also be slow at runtime, as code using these kmers of mixed length will be type unstable.

To solve this, Kmers.jl includes the `DynamicKmer` type.
This type is an immutable, bitstype biosequence, similar to the `Kmer` type, but with the length stored as a run time value.

`DynamicKmer`s are parameterized `DynamicKmer{A <: Alphabet, U <: Unsigned}`, and stores the sequence in a single integer of type `U`. This puts an upper limit on the number elements that can fit in dynamic kmer.
Because the integer also needs to store its runtime length, computing the maximal number of elements in a given
concrete `DynamicKmer` is not straightforward.
It can be obtained with `capacity`.

```@docs
DynamicKmer
capacity
```

For convenience, the aliases `DynamicDNAKmer`, `DynamicRNAKmer` and `DynamicAAKmer` are provided which alias `DynamicKmer{DNAAlphabet{2}, U} where {U <: Unsigned}`, and similar for `RNAAlphabet{2}` and `AminoAcidAlphabet`.

```@docs
DynamicDNAKmer
DynamicRNAKmer
DynamicAAKmer
```

Dynamic kmers can be constructed with the `dmer_str` macro, similar to kmers:

```@docs
@dmer_str
```

### Using dynamically sized kmers
Dynamic kmers behave similar to other biosequences.
Like `LongSequence` and kmers, they can be constructed from strings, bytes and other biosequences.

```jldoctest
julia> m = DynamicRNAKmer{UInt32}("AUGUCGA")
7nt DynamicRNAKmer{UInt32}:
AUGUCGA

julia> complement(m)
7nt DynamicRNAKmer{UInt32}:
UACAGCU

julia> m2 = push(m, DNA_C) # NB: Immutable, so no `push!`
8nt DynamicRNAKmer{UInt32}:
AUGUCGAC
```

### Translating dynamically sized kmers
Dynamic kmers can be translated to obtain `DynamicAAKmers{U}` with various integer types `U`.
The type of `U` is chosen depending on the input type, to ensure that the result will fit in
the output type.

Note that `DynamicDNAKmer{UInt128}`, and its RNA equivalent, contains too many symbols to fit in
`DynamicAAKmer{128}`, and therefore these cannot be translated.

```@docs
BioSequences.translate(::DynamicKmer{<:Union{DNAAlphabet, RNAAlphabet}})
```
