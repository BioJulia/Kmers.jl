module RandomExt

using Random: Random, AbstractRNG, rand, default_rng, Sampler, SamplerTrivial, SamplerType, shuffle

using BioSequences:
    Alphabet,
    @aa_str,
    AminoAcidAlphabet,
    bits_per_symbol,
    iscomplete,
    DefaultAASampler,
    encode,
    symbols

using Kmers:
    Kmer,
    Oligo,
    capacity,
    from_integer,
    get_mask,
    ksize,
    unsafe,
    zero_kmer,
    shift,
    shift_encoding,
    FourBit,
    n_coding_elements,
    derive_type,
    left_shift,
    right_shift,
    _new_dynamic_kmer

# TODO: Add back this AA sampler and the fourbit function to BioSequences.jl
const PROTEOGENIC_AA_ENCODINGS = let
    mem = Vector{UInt8}(undef, 20)
    for (i, sym) in enumerate(aa"ACDEFGHIKLMNPQRSTVWY")
        mem[i] = encode(AminoAcidAlphabet(), sym)
    end
    Tuple(mem)
end

# TODO: Should this be in the core package?
maybe_derive_type(T::Type{Kmer{A, K, N}}) where {A, K, N} = T
maybe_derive_type(T::Type{Kmer{A, K}}) where {A, K} = derive_type(T)

function Random.rand(rng::AbstractRNG, s::SamplerTrivial{T}) where {T <: Kmer}
    iszero(ksize(T)) && throw(ArgumentError("collection must be non-empty"))
    kmer = s[]
    return @inbounds kmer[rand(rng, 1:length(kmer))]
end

function Random.rand(rng::AbstractRNG, s::SamplerTrivial{T}) where {T <: Oligo}
    oligomer = s[]
    isempty(oligomer) && throw(ArgumentError("collection must be non-empty"))
    return @inbounds oligomer[rand(rng, 1:length(oligomer))]
end

"""
    shuffle([rng::AbstractRNG], oligomer::T)::T where {T <: Oligo}

Return an `Oligo` of the same type, length, and symbols as `oligomer`, but
with the symbols randomly permuted by `rng`, which defaults to `Random.default_rng()`.
"""
function Random.shuffle(rng::AbstractRNG, oligomer::Oligo{A, U}) where {A, U}
    n = length(oligomer)
    n < 2 && return oligomer

    bps = bits_per_symbol(A())
    iszero(bps) && return oligomer
    symbol_mask = (one(U) << bps) - one(U)
    u = oligomer.x
    nbits = 8 * sizeof(U)

    # Fisher-Yates shuffle, operating directly on the packed encodings.
    for i in 2:n
        j = rand(rng, 1:i)
        i == j && continue
        ishift = nbits - i * bps
        jshift = nbits - j * bps
        encoding_a = right_shift(u, ishift)
        encoding_b = right_shift(u, jshift)
        # Compute the the xor of the two symbols
        difference = (encoding_a ⊻ encoding_b) & symbol_mask
        # Then xor this difference into the positions of i and j.
        # This switches the symbol at position i to position of symbol j and vice versa
        u = u ⊻ (left_shift(difference, ishift) | left_shift(difference, jshift))
    end
    return _new_dynamic_kmer(A, u)
end

Random.shuffle(oligomer::Oligo) = shuffle(default_rng(), oligomer)

function Random.rand(rng::AbstractRNG, T::Type{<:Oligo{A, U}}, len::Integer) where {A, U}
    # This branch here is crucial to keep, because downstream functions assume no zero-length oligos.
    # These require special handling elsewhere, e.g. in shift_to.
    iszero(len) && return _new_dynamic_kmer(A, zero(U))
    1 <= len <= capacity(T) || throw(ArgumentError("length must be in the 0:capacity(T)"))
    return random_oligomer(rng, T, UInt(len)::UInt)
end

Random.rand(T::Type{<:Oligo}, len::Integer) = rand(default_rng(), T, len)

@inline function random_oligomer(rng::AbstractRNG, T::Type{<:Oligo{A, U}}, len::UInt) where {A <: Alphabet, U}
    return random_oligomer(rng, T, len, iscomplete(Alphabet(T)))
end

@inline function random_oligomer(
        rng::AbstractRNG, T::Type{<:Oligo{A, U}}, len::UInt, ::Val{true}
    ) where {A, U}
    bits = bits_per_symbol(T) * len
    u = shift_to(rand(rng, U), bits) | (len % U)
    return _new_dynamic_kmer(A, u)
end

# Fallback for generic alphabets
function random_oligomer(
        rng::AbstractRNG, T::Type{<:Oligo{A, U}}, len::UInt, ::Val{false}
    ) where {A, U}
    u = zero(U)
    bps = bits_per_symbol(T)
    syms = symbols(A())
    isempty(syms) && throw(ArgumentError("Alphabet cannot be empty"))
    for _ in 1:(len % Int)
        u = left_shift(u, bps) | (encode(A(), rand(rng, syms)) % U)
    end
    return _new_dynamic_kmer(A, shift_to(u, len * bps) | (len % U))
end

@inline function random_oligomer(
        rng::AbstractRNG, T::Type{<:Oligo{A, U}}, len::UInt
    ) where {A <: FourBit, U}
    bits = bits_per_symbol(T) * len
    u = shift_to(random_fourbit_encoding(rng, U), bits) | (len % U)
    return _new_dynamic_kmer(A, u)
end

@inline function random_oligomer(
        rng::AbstractRNG, T::Type{<:Oligo{AminoAcidAlphabet, U}}, len::UInt
    ) where {U}
    u = zero(U)
    bps = bits_per_symbol(T)
    for _ in UInt(1):len
        u |= rand(rng, PROTEOGENIC_AA_ENCODINGS) % U
        u <<= bps
    end
    return _new_dynamic_kmer(AminoAcidAlphabet, shift_to(u, len * bps) | (len % U))
end

function Random.rand(rng::AbstractRNG, ::SamplerType{T}) where {T <: Kmer}
    Tc = maybe_derive_type(T)
    isconcretetype(Tc) || throw(ArgumentError("Cannot sample from abstract K-mer type"))
    return random_kmer(rng, Tc)
end

@inline function random_kmer(rng::AbstractRNG, T::Type{<:Kmer})
    return random_kmer(rng, T, iscomplete(Alphabet(T)))
end

@inline function random_kmer(rng::AbstractRNG, T::Type{<:Kmer{N}}) where {N <: FourBit}
    nce = n_coding_elements(T)
    iszero(nce) && return zero_kmer(T)
    tail = ntuple(i -> random_fourbit_encoding(rng, UInt), nce - 1)
    head = random_fourbit_encoding(rng, UInt) & get_mask(T)
    return T(unsafe, (head, tail...))
end

@inline function random_kmer(rng::AbstractRNG, T::Type{<:Kmer{AminoAcidAlphabet}})
    kmer = zero_kmer(T)
    for _ in 1:ksize(T)
        kmer = shift_encoding(kmer, rand(rng, PROTEOGENIC_AA_ENCODINGS) % UInt)
    end
    return kmer
end

@inline function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, ::Val{true})
    bits = bits_per_symbol(T) * ksize(T)
    return T(unsafe, random_tuples(rng, Val(bits)))
end

function random_kmer(rng::AbstractRNG, T::Type{<:Kmer}, ::Val{false})
    letters = symbols(Alphabet(T))
    isempty(letters) && throw(ArgumentError("Alphabet cannot be empty"))
    kmer = zero_kmer(T)
    for _ in 1:ksize(T)
        kmer = shift(kmer, rand(rng, letters))
    end
    return kmer
end

@inline function random_tuples(rng::AbstractRNG, v::Val{bits}) where {bits}
    usize = 8 * sizeof(UInt)
    return random_tuples(rng, v, Val{iszero(mod(bits, usize))}())
end

@inline function random_tuples(rng::AbstractRNG, ::Val{bits}, ::Val{true}) where {bits}
    return ntuple(i -> rand(rng, UInt), cld(bits, 8 * sizeof(UInt)))
end

@inline function random_tuples(rng::AbstractRNG, ::Val{bits}, ::Val{false}) where {bits}
    usize = 8 * sizeof(UInt)
    tail = ntuple(i -> rand(rng, UInt), cld(bits, usize) - 1)
    head = rand(rng, UInt) & (UInt(1) << mod(bits, usize) - UInt(1))
    return (head, tail...)
end

function random_fourbit_encoding(rng::AbstractRNG, ::Type{U}) where {U <: Unsigned}
    mask = rand(rng, U)
    enc = div(typemax(U), 0x0f)
    enc += enc & mask
    return enc + 0x03 * (enc & (mask >>> 1))
end

shift_to(u::Unsigned, bits::UInt) = left_shift(u, (8 * sizeof(u)) % UInt - bits)

end # module
