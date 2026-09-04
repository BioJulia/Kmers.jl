```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```

## Dynamically sized kmers

### Quick Start

```jldoctest
julia> d = dmer"TAGCAT"d  # Create DNA kmer from string literal
6nt DNAOligo{UInt64}:
TAGCAT

julia> reverse_complement(d)
6nt DNAOligo{UInt64}:
ATGCTA

julia> push(d, DNA_G)  # Returns new instance (immutable)
7nt DNAOligo{UInt64}:
TAGCATG
```

### Overview

Sometimes, one requires using kmers of differing lengths in the same workload.
An example could be representing primers, which can be of length 18-24.
Here, using the `Kmer` type would cause code to specialize on each kmer length.
Besides causing both excessive compilation and code generation, it will also be slow at runtime, as code using these kmers of mixed length will be type unstable.

To solve this, Kmers.jl includes the `Oligo` type.
This type is an immutable, bitstype biosequence, similar to the `Kmer` type, but with the length stored as a runtime value rather than a compile-time type parameter.

```@docs
Oligo
```

### Basic Properties

`Oligo` has several important characteristics:

* **Immutable**: All operations return new instances. There is no `push!`, only `push`.
* **Bitstype**: Stored inline in a single unsigned integer.
* **Runtime length**: Unlike `Kmer`, the length is not part of the type, avoiding type instability for variable-length workloads.
* **Size limits**: Each `Oligo` type has a maximum capacity determined by its alphabet and backing integer type.
* **Performance**: Slightly slower than `Kmer` but much faster than `LongSequence` for small sequences.

### Type Parameters

`Oligo{A, U}` is parameterized by `A`, its `Alphabet`, and `U`, the backing unsigned integer type.
Thus, an `Oligo{DNAAlphabet{4}, UInt32}` is 4 bytes in size and contains 4-bit DNA.

For interactive use, `Oligo{A}(source)` and the `DNAOligo(source)`,
`RNAOligo(source)`, and `AAOligo(source)` aliases use a fixed `UInt64`
backing type. This keeps the result type independent of the source's runtime
length.

The backing integer `U` stores both the sequence data and the runtime length.
This imposes capacity limits. Use `capacity(T)` to determine the maximum number of symbols for a given type.

```@docs
capacity
```

For convenience, type aliases are provided:

```@docs
DNAOligo
RNAOligo
AAOligo
```

### Capacity and Size Limits

The maximum number of symbols depends on both the alphabet and the backing integer type:

| Type | Capacity | Bits per symbol | Notes |
|------|----------|-----------------|-------|
| `DNAOligo{UInt32}` | 14 | 2 | Good for short primers |
| `DNAOligo{UInt64}` | 29 | 2 | Standard choice for DNA/RNA |
| `RNAOligo{UInt64}` | 29 | 2 | Same capacity as DNA |
| `AAOligo{UInt128}` | 15 | 8 | Limited by bytes per symbol |

Choose a larger backing integer for longer sequences, but be aware that integers larger than 64 bits
typically become slower as they get larger.
For very large dynamic kmers, you can use the `BitIntegers.jl` package that provides, e.g., `UInt512`,
but make sure to test whether, for your application, huge dynamic kmers are still faster than
`BioSequences.LongSequence`.

### Construction Methods

Dynamic kmers can be constructed with the `@dmer_str` macro, similar to kmers:

```@docs
@dmer_str
```

Like other `BioSequence` types, they can also be constructed from strings, `AbstractVector{UInt8}`
(interpreted as containing ASCII), and other `BioSequence`s.

```jldoctest
julia> DNAOligo("TAGCAT")  # From string; uses UInt64
6nt DNAOligo{UInt64}:
TAGCAT

julia> RNAOligo(rna"AUGCUA")  # From BioSequence
6nt RNAOligo{UInt64}:
AUGCUA

julia> DNAOligo([DNA_T, DNA_A, DNA_G])  # From iterable
3nt DNAOligo{UInt64}:
TAG

julia> AAOligo([0x61, 0x63])  # From ASCII AbstractVector{UInt8}
2aa AAOligo{UInt64}:
AC
```

For extracting a runtime-length window without bounds checks, use the public
five-argument [`unsafe_extract`](@ref) overload. Unlike the `Kmer` overload, its
runtime `len` argument specifies the result length. The source range and
destination capacity must be checked by the caller; symbol validity is still
checked during recoding.

```jldoctest
julia> source = codeunits("CCTAGCAA");

julia> Kmers.unsafe_extract(Kmers.AsciiEncode(), DNAOligo{UInt32}, source, 3, 4)
4nt DNAOligo{UInt32}:
TAGC
```

### Common Operations

Dynamic kmers behave similarly to other biosequences, supporting biological transformations, indexing, and length-changing operations.

#### Biological Transformations

```jldoctest
julia> d = dmer"TAGCAT"d
6nt DNAOligo{UInt64}:
TAGCAT

julia> reverse_complement(d)
6nt DNAOligo{UInt64}:
ATGCTA

julia> canonical(d)
6nt DNAOligo{UInt64}:
ATGCTA
```

#### Modifying Length

All operations return new instances since `Oligo` is immutable.
They use, e.g., `pop` instead of `pop!`.
Instead of `popfirst!` and `pushfirst!`, use `pop_first` and `push_first`
(note the underscore):

```jldoctest
julia> d = dmer"TAG"d
3nt DNAOligo{UInt64}:
TAG

julia> push(d, DNA_C)  # Add to end
4nt DNAOligo{UInt64}:
TAGC

julia> push_first(d, DNA_C)  # Add to beginning
4nt DNAOligo{UInt64}:
CTAG

julia> shift(d, RNA_U)  # Append while discarding the first symbol
3nt DNAOligo{UInt64}:
AGT

julia> shift_first(d, RNA_U)  # Prepend while discarding the last symbol
3nt DNAOligo{UInt64}:
TTA
```

#### Indexing and Slicing
Slicing returns a value of the same type. Unlike with `Kmer`, this is type stable.
Use `Base.setindex` to create a new kmer with a given `BioSymbol` replaced.

```jldoctest
julia> d = dmer"TAGCAT"d
6nt DNAOligo{UInt64}:
TAGCAT

julia> d[2:4]
3nt DNAOligo{UInt64}:
AGC

julia> d[[6, 2, 2, 1]]
4nt DNAOligo{UInt64}:
TAAT

julia> d[[true, false, true, false, false, true]]
3nt DNAOligo{UInt64}:
TGT

julia> Base.setindex(d, 'G', 2)  # Immutable, returns new kmer
6nt DNAOligo{UInt64}:
TGGCAT
```

#### Integer Conversion
Like `Kmer`s, `Oligo` can be converted to and from integers.
Unlike the `Kmer` method, length is required when using `from_integer`. The input
may have a different unsigned integer width from the result's backing type; only
the low coding bits are used:

```@docs
as_integer(::Oligo)
from_integer(T::Type{Oligo{A, U}}, x::Unsigned, len::Int) where {A <: Alphabet, U <: Unsigned}
```

```jldoctest
julia> original = DNAOligo{UInt32}("TAGC");

julia> encoded = UInt128(as_integer(original));

julia> from_integer(DNAOligo{UInt64}, encoded, length(original))
4nt DNAOligo{UInt64}:
TAGC
```

### Type Conversions and Compatibility

#### Between Backing Integer Types

Dynamic kmers can be converted between different backing integer types. Comparisons
require matching backing types, so widen the smaller value before comparing:

```jldoctest
julia> d32 = DNAOligo{UInt32}("TAGC")
4nt DNAOligo{UInt32}:
TAGC

julia> d64 = DNAOligo{UInt64}(d32)  # Widen to larger type
4nt DNAOligo{UInt64}:
TAGC

julia> d64 == DNAOligo{UInt64}("TAGC")
true
```

#### Comparison Compatibility

Equality, ordering, and `cmp` are defined for `Oligo`s with the same backing
integer type and either the same alphabet or compatible DNA/RNA alphabets of the
same bit width. Compatible DNA and RNA values compare by their shared encoding and
equal values have equal hashes:

```jldoctest
julia> dna = DNAOligo{UInt32}("TAGC"); rna = RNAOligo{UInt32}("UAGC");

julia> dna == rna
true

julia> hash(dna) == hash(rna)
true
```

Comparing different backing widths, different bits-per-symbol encodings, or an
`Oligo` with another `BioSequence` type throws a `MethodError`. Convert both
values to a common `Oligo` type first.

#### Translating Dynamic Kmers

Dynamic kmers can be translated to obtain `AAOligo{U}` with various integer types `U`.
The type of `U` is chosen depending on the input type to ensure that the result will fit in
the output type.

By default, the largest `AAOligo` type is `AAOligo{UInt128}`. However, if the package
BitIntegers.jl is loaded, Kmers.jl will make use of larger integer sizes,
currently up to `UInt1024`.

```@docs
BioSequences.translate(::Oligo{<:Union{DNAAlphabet, RNAAlphabet}})
```

### Intentional Differences from `Kmer`

The following `Kmer` facilities are intentionally not generalized to `Oligo`:

* Kmers.jl does not provide `Oligo`-specific sequence iterators. The existing
  iterator types produce fixed-length `Kmer`s. However, construction of `Oligo`
  from `Kmer` is cheap, so to construct sliding window `Oligo`, you can do
  something like `(DNAOligo{UInt32}(i) for i in FwDNAMers{3}("TAGTGCA"))`.
