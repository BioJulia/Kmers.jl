module RandomExt

using Random: Random, AbstractRNG, rand, default_rng, Sampler, SamplerTrivial, SamplerType

using BioSequences: Alphabet,
    @aa_str,
    AminoAcidAlphabet,
    bits_per_symbol,
    iscomplete,
    DefaultAASampler,
    encode

using Kmers: Kmer,
    get_mask,
    ksize,
    unsafe,
    zero_kmer,
    shift,
    shift_encoding,
    FourBit,
    n_coding_elements,
    derive_type

# TODO: Add back this AA sampler and the fourbit function to BioSequences.jl
const PROTEOGENIC_AA_ENCODINGS = let
    mem = Memory{UInt8}(undef, 20)
    for (i, sym) in enumerate(aa"ACDEFGHIKLMNPQRSTVWY")
        mem[i] = encode(AminoAcidAlphabet(), sym)
    end
    mem
end

# TODO: Should this be in the core package?
maybe_derive_type(T::Type{Kmer{A, K, N}}) where {A, K, N} = T
maybe_derive_type(T::Type{Kmer{A, K}}) where {A, K} = derive_type(T)

function Random.rand(rng::AbstractRNG, s::SamplerTrivial{T}) where {T <: Kmer}
    iszero(ksize(T)) && throw(ArgumentError("collection must be non-empty"))
    kmer = s[]
    return @inbounds kmer[rand(rng, 1:length(kmer))]
end

function Random.rand(rng::AbstractRNG, ::SamplerType{T}) where {T <: Kmer}
    Tc = maybe_derive_type(T)
    isconcretetype(Tc) || throw(ArgumentError("Cannot sample from abstract K-mer type"))
    return random_kmer(rng, Tc)
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer})
    return random_kmer(rng, T, iscomplete(Alphabet(T)))
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer{N}}) where {N <: FourBit}
    nce = n_coding_elements(T)
    iszero(nce) && return zero_kmer(T)
    tail = ntuple(i -> random_fourbit_encoding(), nce - 1)
    head = random_fourbit_encoding() & get_mask(T)
    return T(unsafe, (head, tail...))
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer{AminoAcidAlphabet}})
    kmer = zero_kmer(T)
    for _ in 1:ksize(T)
        kmer = shift_encoding(kmer, rand(PROTEOGENIC_AA_ENCODINGS) % UInt)
    end
    return kmer
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, ::Val{true})
    bits = bits_per_symbol(T) * ksize(T)
    return T(unsafe, random_tuples(Val(bits)))
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, ::Val{false})
    letters = symbols(Alphabet(T))
    isempty(letters) && throw(ArgumentError("Alphabet cannot be empty"))
    kmer = zero_kmer(T)
    for i in 1:ksize(T)
        kmer = shift(kmer, rand(rng, letters))
    end
    return kmer
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

end # module
