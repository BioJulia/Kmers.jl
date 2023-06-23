```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```

!!! warning
    The value of hashes are guaranteed to be reproducible for a given version
    of Kmers.jl and Julia, but may __change__ in new minor versions of Julia
    or Kmers.jl

## Hashing
Kmers.jl implements `Base.hash`, yielding a `UInt` value:

```jldoctest; filter = r"^0x[0-9a-fA-F]+$"
julia> hash(mer"UGCUGUAC"r)
0xe5057d38c8907b22
```

The implementation of `Base.hash` for kmers strikes a compromise between providing a high-quality (non-cryptographic) hash, while being reasonably fast.
While hash collisions can easily be found, they are unlikely to occur at random.
When kmers are of the same (or compatible) alphabets, different kmers hash to different values
(not counting the occational hash collision), even when they have the same underlying bitpattern:

```jldoctest
julia> using BioSequences: encoded_data

julia> a = mer"TAG"d; b = mer"AAAAAAATAG"d;

julia> encoded_data(a) === encoded_data(b)
true

julia> hash(a) == hash(b)
false
```

When they are of compatible alphabets, and have the same content, they hash to the same value.
Currently, only DNA and RNA of the alphabets `DNAAlphabet` and `RNAAlphabet` are compatible:

```jldoctest
julia> a = mer"UUGU"r; b = mer"TTGT"d;

julia> a == b # equal
true

julia> a === b # not egal
false

julia> hash(a) === hash(b)
true
```

For some applications, fast hashing is absolutely crucial. For these cases, Kmers.jl provides [`fx_hash`](@ref), which trades off hash quality for speed:

```@docs
fx_hash
```