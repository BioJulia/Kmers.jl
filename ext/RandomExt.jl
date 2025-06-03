module RandomExt

using Random: Random, AbstractRNG, rand, default_rng
using BioSequences: bits_per_symbol, iscomplete, SamplerUniform, Alphabet
using Kmers: Kmer,
    get_mask,
    ksize,
    unsafe,
    zero_kmer,
    shift,
    FourBit,
    n_coding_elements,
    derive_type

maybe_derive_type(T::Type{Kmer{A, K, N}}) where {A, K, N} = T
maybe_derive_type(T::Type{Kmer{A, K}}) where {A, K} = derive_type(T)

Random.rand(T::Type{<:Kmer}) = rand(default_rng(), T)

function Random.rand(rng::AbstractRNG, T::Type{<:Kmer{<:FourBit}})
    T = maybe_derive_type(T)
    nce = n_coding_elements(T)
    iszero(nce) && return zero_kmer(T)
    tail = ntuple(i -> random_fourbit_encoding(), nce - 1)
    head = random_fourbit_encoding() & get_mask(T)
    return T(unsafe, (head, tail...))
end

function Random.rand(rng::AbstractRNG, T::Type{<:Kmer})
    return random_kmer(rng, T, iscomplete(Alphabet(T)))
end

function random_tuples(::Val{bits}) where {bits}
    iszero(bits) && return ()
    usize = 8 * sizeof(UInt)
    tail = ntuple(i -> rand(UInt), cld(bits, usize) - 1)
    head = rand(UInt) & (UInt(1) << mod(bits, usize) - UInt(1))
    return (head, tail...)
end

function random_fourbit_encoding()
    enc = 0x1111111111111111 % UInt
    mask = rand(UInt)
    enc = ((enc & mask) << 1) | (enc & ~mask)
    mask >>>= 1
    return ((enc & mask) << 2) | (enc & ~mask)
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, ::Val{true})
    bps = bits_per_symbol(T) * ksize(T)
    T = maybe_derive_type(T)
    return T(unsafe, random_tuples(Val(bps)))
end

function random_kmer(rng::AbstractRNG, AbstractRNG, T::Type{<:Kmer}, ::Val{false})
    A = Alphabet(T)
    letters = symbols(A)
    sampler = SamplerUniform{eltype(A)}(letters)
    return random_kmer(rng, T, sampler)
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, sampler::SamplerUniform)
    kmer = zero_kmer(T)
    for i in 1:ksize(T)
        kmer = shift(kmer, rand(rng, sampler))
    end
    return kmer
end


end # module
