```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```

## The `Kmer` type
The central type of Kmers.jl is the `Kmer`.
A `Kmer` is an immutable, bitstype `BioSequence`, with a length known at compile
time. Compared to `LongSequence` in BioSequences.jl,
this gives to one advantage, and comes with two disadvantages:
* Kmers are much faster than `LongSequence`, as they can be stored in registers.
* As kmers gets longer, the code gets increasingly inefficient, as the unrolling
  and inlining of the immutable operations breaks down.
* Since their length is part of their type, any operation that results in a kmer
  whose length cannot be determined at compile time will be type unstable.
  This includes slicing a kmer, pushing and popping it, and other operations.

The `Kmer` type is (roughly) defined as
```julia
struct Kmer{A <: Alphabet, K, N} <: BioSequence{A}
    x::NTuple{N, UInt}
end
```
Where:
* `A` is the `Alphabet` as defined in BioSequences.jl.
* `K` is the length.
* `N` is an extra type parameter derived from the first two,
  which exists only because Julia does not allow computed type parameters.

### Construction
Kmers can be constructed from a `BioSequence` or `AbstractString` by explicitly
specifying the length of the sequence:

```jldoctest
julia> Kmer{DNAAlphabet{2}, 5, 1}("TAGCT")
DNA 5-mer:
TAGCT
```

The final type parameter can be elided, in which case it will be inferred:

```jldoctest
julia> Kmer{DNAAlphabet{2}, 5}("TAGCT")
DNA 5-mer:
TAGCT
```

Kmers with alphabets `DNAAlphabet{2}`, `RNAAlphabet{2}` and `AminoAcidAlphabet`
can be created with the type aliases `DNAKmer`, `RNAKmer` and `AAKmer`:

```jldoctest
julia> DNAKmer{3}("tag")
DNA 3-mer:
TAG

julia> AAKmer{5}("PWYSK")
AminoAcid 5-mer:
PWYSK
```

For kmers with an `Alphabet` that implement `BioSequences.AsciiAlphabet`, they can also be constructed from `AbstractVector{UInt8}`, in which case the vector is interpreted as being bytes of ASCII text:

```jldoctest
julia> AAKmer{3}([0x65, 0x67, 0x7a])
AminoAcid 3-mer:
EGZ
```

When constructing from an `AbstractString` (or byte vector), uracil (`U`) and thymine `T` are treated differently - a `U` cannot be read as thymine:

```jldoctest
julia> DNAKmer{3}("UAG")
ERROR: cannot encode 0x55 (Char 'U') in DNAAlphabet{2}
[...]
```

However, when constructing from a `BioSequence`, these nucleotides are considered
interchangeable:

```jldoctest
julia> RNAKmer{4}(dna"TATC")
RNA 4-mer:
UAUC
```

Finally, kmers can be constructed with a string literal `@mer_str`, where the string must be appended with `d` for DNA, `r` for RNA, or `a` for amino acid:

```jldoctest
julia> mer"UGCUGA"r
RNA 6-mer:
UGCUGA

julia> mer"EDEHL"a
AminoAcid 5-mer:
EDEHL
```

Since the literals produce the kmer at parse time and inserts it directly into the parsed code, this will always be type stable,
and the overhead related to parsing the string will not be paid:

```jldoctest; filter = [r"^\s*0\.\d+ seconds.+"s, r"^\d+"s]
julia> function count_aaas(dna)
           x = 0
           for kmer in FwDNAMers{3}(dna)
               # The parsing happens once here, when the
               # code is parsed, and is fine to have in the loop
               x += kmer == mer"AAA"d
           end
           x
       end;

julia> seq = randseq(DNAAlphabet{2}(), 100_000_000);

julia> @time count_aaas(seq)
  0.193463 seconds (32.05 k allocations: 2.051 MiB, 21.88% compilation time)
1563330
```


### Indexing
Kmers.jl supports most normal indexing, such as scalar indexing:

```jldoctest
julia> mer"CAGCU"r[3]
RNA_G
```

Slicing

```jldoctest
julia> mer"AGGCTA"d[2:5]
DNA 4-mer:
GGCT
```

And indexing with boolean vectors, and vectors of indices:

```jldoctest
julia> m = mer"MDGKRY"a;

julia> m[[true, false, true, true, false, true]]
AminoAcid 4-mer:
MGKY

julia> m[[4,2]]
AminoAcid 2-mer:
KD
```

### A note on type stability
!!! warning
    Except scalar indexing which always returns a single symbol, all the operations
    above are _type unstable_, since the length (and thus type) of the resulting 
    kmer depend on the input value, not its type.

However, type unstable functions may be type-stable, if the indexing value is
known at compile time, and the Julia compiler uses constant folding:

```jldoctest
julia> f(x) = x[2:5]; # 2:5 is a compile time constant

julia> Test.@inferred f(mer"UCGUAGC"r)
RNA 4-mer:
CGUA
```

### Reference
```@docs
Kmer
Mer
@mer_str
DNAKmer
RNAKmer
AAKmer
DNACodon
RNACodon
pop
pop_first
push
push_first
shift
shift_first
```