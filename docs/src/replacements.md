```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
    using Random
end
```

# [Building kmer replacements](@id replacements)
_Kmer replacements_ is the general term for sequences that can be represented
computationally like kmers, but are sampled differently. Examples include minimizers, strobemers, syncmers and k-min-mers.

Since there is no end to the variations of kmer replacements, Kmers.jl does not try to implement all of them.
Instead, Kmers.jl implements the base kmer type, and exposes some efficient primitives to allow downstream users to build kmer replacements.

These functions are:

```@docs
unsafe_extract
shift_encoding
unsafe_shift_from
```

# Example: Minimizers
Minimizers are currently the most common kmer replacement.
They are defined as the minimum of W consecutive kmers, as ordered by some ordering O.

If we use [`fx_hash`](@ref) as the ordering function, and assume K and W are known at compile time, we can implement it reasonably efficiently like so:

```jldoctest
function unsafe_extract_minimizer(
    seq::LongDNA{2},
    i::Int,
    ::Val{K},
    ::Val{W},
) where {K, W}
    T = derive_type(Kmer{DNAAlphabet{2}, K})
    kmer = Kmers.unsafe_extract(Kmers.Copyable(), T, seq, i)
    hash = fx_hash(kmer)
    for offset in 0:W-2
        new_kmer = Kmers.unsafe_shift_from(Kmers.Copyable(), kmer, seq, i+K+offset, Val(1))
        new_hash = fx_hash(new_kmer)
        if new_hash < hash
            hash = new_hash
            kmer = new_kmer
        end
    end
    kmer
end

rng = Random.Xoshiro(1)
unsafe_extract_minimizer(randseq(rng, DNAAlphabet{2}(), 100), 1, Val(5), Val(9))

# output
DNA 5-mer:
TATCA
```
