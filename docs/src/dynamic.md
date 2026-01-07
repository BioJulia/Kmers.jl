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

`DynamicKmer`s are parameterized `DynamicKmer{A <: Alphabet, U <: Unsigned}`, and stores the sequence in a single integer of type `U`. This puts an upper limit on the number of coding bits.
For example, here is the maximum number of coding bits when using `A = DNAAlphabet{2}`:

| U       | Available bits  |
|:--------|:----------------|
| UInt8   | 6               |
| UInt16  | 12              |
| UInt32  | 28              |
| UInt64  | 58              |
| UInt128 | 122             |

The remaining bits are used to store the length of the kmer.

Dynamic kmers can be constructed with the `dmer_str` macro, similar to kmers:

```@docs
@dmer_str
```

### Using dynamically sized kmers
* Dynamic kmers can be constructed from a normal BioSequence or string, and can then be treated like a normal `BioSequence`.

```jldoctest
julia> m = DynamicRNAKmer{UInt32}("AUGUCGA")
7nt DynamicRNAKmer{UInt32}:
AUGUCGA

julia> complement(m)
7nt DynamicRNAKmer{UInt32}:
UACAGCU
```

```@docs
DynamicKmer
DynamicDNAKmer
DynamicRNAKmer
DynamicAAKmer
```