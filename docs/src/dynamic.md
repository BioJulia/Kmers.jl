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
6nt DNAOligomer{UInt64}:
TAGCAT

julia> reverse_complement(d)
6nt DNAOligomer{UInt64}:
ATGCTA

julia> push(d, DNA_G)  # Returns new instance (immutable)
7nt DNAOligomer{UInt64}:
TAGCATG
```

### Overview

Sometimes, one requires using kmers of differing lengths in the same workload.
An example could be representing primers, which can be of length 18-24.
Here, using the `Kmer` type would cause code to specialize on each kmer length.
Besides causing both excessive compilation and code generation, it will also be slow at runtime, as code using these kmers of mixed length will be type unstable.

To solve this, Kmers.jl includes the `Oligomer` type.
This type is an immutable, bitstype biosequence, similar to the `Kmer` type, but with the length stored as a runtime value rather than a compile-time type parameter.

```@docs
Oligo
```

### Basic Properties

`Oligomer` has several important characteristics:

- **Immutable**: All operations return new instances. There is no `push!`, only `push`.
- **Bitstype**: Stored inline in a single unsigned integer.
- **Runtime length**: Unlike `Kmer`, the length is not part of the type, avoiding type instability for variable-length workloads.
- **Size limits**: Each `Oligomer` type has a maximum capacity determined by its alphabet and backing integer type.
- **Performance**: Slightly slower than `Kmer` but much faster than `LongSequence` for small sequences.

### Type Parameters

`Oligomer{A, U}` is parameterized by `A`, its `Alphabet`, and `U`, the backing unsigned integer type.
Thus, a `Oligomer{DNAAlphabet{4}, UInt32}` is 4 bytes in size, and contains 4-bit DNA.

The backing integer `U` stores both the sequence data and the runtime length.
This imposes capacity limits. Use `capacity(T)` to determine the maximum number of symbols for a given type.

```@docs
capacity
```

For convenience, type aliases are provided:

```@docs
DNAOligomer
RNAOligomer
AAOligomer
```

### Capacity and Size Limits

The maximum number of symbols depends on both the alphabet and the backing integer type:

| Type | Capacity | Bits per symbol | Notes |
|------|----------|-----------------|-------|
| `DNAOligomer{UInt32}` | 14 | 2 | Good for short primers |
| `DNAOligomer{UInt64}` | 30 | 2 | Standard choice for DNA/RNA |
| `RNAOligomer{UInt64}` | 30 | 2 | Same capacity as DNA |
| `AAOligomer{UInt128}` | 15 | 8 | Limited by byte-per-symbol |

Choose a larger backing integer for longer sequences, but be aware that integers larger than 64 bits
typically become slower the larger they are.
For very large dynamic kmers, you can use the `BitIntegers.jl` package that provide e.g. `UInt512`,
but make sure to test that for your application, huge dynamic kmers are still faster than
`BioSequences.LongSequence`. 

### Construction Methods

Dynamic kmers can be constructed with the `@dmer_str` macro, similar to kmers:

```@docs
@dmer_str
```

Like other `BioSequence`s, they can also be constructed from strings, `AbstractVector{UInt8}`
(interpreted as containing ASCII), and other `BioSequence`s.

```jldoctest
julia> DNAOligomer{UInt64}("TAGCAT")  # From string
6nt DNAOligomer{UInt64}:
TAGCAT

julia> RNAOligomer{UInt64}(rna"AUGCUA")  # From BioSequence
6nt RNAOligomer{UInt64}:
AUGCUA

julia> DNAOligomer{UInt64}([DNA_T, DNA_A, DNA_G])  # From iterable
3nt DNAOligomer{UInt64}:
TAG

julia> AAOligomer{UInt64}([0x61, 0x63])  # From ASCII AbstractVector{UInt8}
2aa AAOligomer{UInt64}:
AC
```

### Common Operations

Dynamic kmers behave similar to other biosequences, supporting biological transformations, indexing, and length-changing operations.

#### Biological Transformations

```jldoctest
julia> d = dmer"TAGCAT"d
6nt DNAOligomer{UInt64}:
TAGCAT

julia> reverse_complement(d)
6nt DNAOligomer{UInt64}:
ATGCTA

julia> canonical(d)
6nt DNAOligomer{UInt64}:
ATGCTA
```

#### Modifying Length

All operations return new instances since `Oligomer` is immutable.
They use e.g. `pop` instead of `pop!`.
Instead of `popfirst!` and `pushfirst!`, it uses `pop_first` and `push_first`
(note the underscore):

```jldoctest
julia> d = dmer"TAG"d
3nt DNAOligomer{UInt64}:
TAG

julia> push(d, DNA_C)  # Add to end
4nt DNAOligomer{UInt64}:
TAGC

julia> push_first(d, DNA_C)  # Add to beginning
4nt DNAOligomer{UInt64}:
CTAG
```

#### Indexing and Slicing
Slicing returns a value of the same type. Unlike for `Kmer`, this is type stable.
Use Base.setindex to create a new kmer with a given `BioSymbol` replaced.


```jldoctest
julia> d = dmer"TAGCAT"d
6nt DNAOligomer{UInt64}:
TAGCAT

julia> d[2:4]
3nt DNAOligomer{UInt64}:
AGC

julia> Base.setindex(d, 'G', 2)  # Immutable, returns new kmer
6nt DNAOligomer{UInt64}:
TGGCAT
```

#### Integer Conversion
Like `Kmer`s, `Oligomer` can be converted to and from integers.
Unlike the `Kmer` method, length is required when using `from_integer`:

```@docs
as_integer(::Oligomer)
from_integer(T::Type{Oligomer{A, U}}, x::U, len::Int) where {A <: Alphabet, U <: Unsigned}
```

### Type Conversions and Compatibility

#### Between Backing Integer Types

Dynamic kmers can be converted between different backing integer types:

```jldoctest
julia> d32 = DNAOligomer{UInt32}("TAGC")
4nt DNAOligomer{UInt32}:
TAGC

julia> d64 = DNAOligomer{UInt64}(d32)  # Widen to larger type
4nt DNAOligomer{UInt64}:
TAGC

julia> d64 == d32  # Comparable across backing types
true
```

#### Translating Dynamic Kmers

Dynamic kmers can be translated to obtain `AAOligomer{U}` with various integer types `U`.
The type of `U` is chosen depending on the input type, to ensure that the result will fit in
the output type.

By default, the largest `AAOligomer` type is `AAOligomer{UInt128}`. However, if the package
BitIntegers.jl is loaded, Kmers.jl will make use of larger integer sizes,
currently up to `UInt1024`.

```@docs
BioSequences.translate(::Oligomer{<:Union{DNAAlphabet, RNAAlphabet}})
```
