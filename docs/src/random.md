```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using Kmers
end
```

# Generating random sequences

You can generate random kmers using `Base.rand` function.

```@docs
Base.rand(::Type{<:Kmer})
```