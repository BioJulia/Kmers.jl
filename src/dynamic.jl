struct DynamicKmer{A <: Alphabet, U <: Unsigned} <: BioSequence{A}
    # Lower L bits: Length
    # Upper bits, from top to bottom: bits
    # E.g. A = 2bit and U = UInt8, TG is stored:
    #  11 10 00                    10
    #  T  G  unused (always zero)  length

    x::U

    global function _new_dynamic_kmer(::Type{A}, x::U) where {A, U}
        return new{A, U}(x)
    end
end

const DynamicDNAKmer{U} = DynamicKmer{DNAAlphabet{2}, U}
const DynamicRNAKmer{U} = DynamicKmer{RNAAlphabet{2}, U}
const DynamicAAKmer{U} = DynamicKmer{AminoAcidAlphabet, U}

Base.@constprop :aggressive Base.@assume_effects :foldable function max_coding_bits(
        ::Type{DynamicKmer{A, U}}
    ) where {A, U}
    bps = BioSequences.bits_per_symbol(A())
    iszero(bps) && return 0
    B = 8 * sizeof(U)
    for S in div(B, bps):-1:0
        Lb = 64 - leading_zeros(S)
        Sb = bps * S
        Lb + Sb ≤ B && return Sb
    end
    0
end

Base.@constprop :aggressive Base.@assume_effects :foldable function length_bits(
        T::Type{DynamicKmer{A, U}}
    ) where {A, U}
    8 * sizeof(U) - max_coding_bits(T)
end

@inline function length_mask(T::Type{<:DynamicKmer})
    U = BioSequences.encoded_data_eltype(T)
    return one(U) << length_bits(T) - one(U)
end

@inline function top_mask(::Type{U}, len::Integer) where {U <: Unsigned}
    return left_shift(typemax(U), (8 * sizeof(U) - len))
end

@inline function top_mask(T::Type{<:DynamicKmer}, len::Integer)
    return top_mask(BioSequences.encoded_data_eltype(T), len)
end

@inline function coding_bits(x::DynamicKmer)
    return BioSequences.bits_per_symbol(x) * length(x)
end

@inline function noncoding_bits(x::DynamicKmer)
    return 8 * sizeof(x) - coding_bits(x)
end

@inline function coding_mask(x::DynamicKmer)
    return top_mask(typeof(x), coding_bits(x))
end

capacity(T::Type{<:DynamicKmer}) = div(max_coding_bits(T), BioSequences.bits_per_symbol(Alphabet(T)))

BioSequences.encoded_data_eltype(::Type{DynamicKmer{A, U}}) where {A, U} = U

function BioSequences.extract_encoded_element(x::DynamicKmer, i::Integer)
    bps = BioSequences.bits_per_symbol(x)
    shift = 8 * sizeof(x) - (i * bps)
    u = right_shift(x.x, shift)
    mask = one(x.x) << (bps) - one(x.x)
    return u & mask
end

Base.length(x::DynamicKmer) = (x.x & length_mask(typeof(x))) % Int

function DynamicKmer{A, U}(itr) where {A <: Alphabet, U <: Unsigned}
    T = DynamicKmer{A, U}
    bps = BioSequences.bits_per_symbol(T)
    shift = 8 * sizeof(T)
    u = zero(U)
    cap = capacity(T)
    len = 0
    for i in itr
        len += 1
        shift -= bps
        len > cap && error("Iterator size exceeds maximum capacity of dynamic kmer")
        enc = BioSequences.encode(A(), i) % U
        u |= left_shift(enc, shift)
    end
    return _new_dynamic_kmer(A, (len % U) | u)
end

function DynamicKmer{A, U}(kmer::Kmer{A}) where {A <: Alphabet, U <: Unsigned}
    T = DynamicKmer{A, U}
    K = length(kmer)
    if capacity(T) < K
        error("Kmer size exceeds maximum size of dynamic kmer")
    end
    u = as_integer(kmer) % U
    bps = BioSequences.bits_per_symbol(A())
    shift = 8 * sizeof(T) - (K * bps)
    return _new_dynamic_kmer(A, left_shift(u, shift) | (K % U))
end

function Kmer{A, K}(x::DynamicKmer{A}) where {A <: Alphabet, K}
    return derive_type(Kmer{A, K})(x)
end

@assert UInt == UInt64

function Kmer{A, K, N}(x::DynamicKmer) where {A <: Alphabet, K, N}
    length(x) == K || error("Must construct kmer from length K DynamicKmer")
    # This is now a compile time constant
    noncoding = 8 * sizeof(x) - BioSequences.bits_per_symbol(x) * K
    return if N == 0
        Kmer{A, K, N}(unsafe, ())
    elseif N == 1
        u = right_shift(x.x % UInt64, noncoding)
        Kmer{A, K, N}(unsafe, (u,))
    elseif N == 2
        u = right_shift(x.x, noncoding)
        t1 = (u >> 64) % UInt64
        t2 = u % UInt64
        Kmer{A, K, N}(unsafe, (t1, t2))
    else
        error("Unreachable")
    end
end

const HASH_MASK = 0x6ff6e9f0462d5162 % UInt

Base.hash(x::DynamicKmer, h::UInt64) = hash(x.x, h ⊻ HASH_MASK)
Base.:(==)(a::DynamicKmer, b::DynamicKmer) = a.x == b.x
Base.isless(a::DynamicKmer{A}, b::DynamicKmer{A}) where {A} = isless(a.x, b.x)
Base.cmp(a::DynamicKmer{A}, b::DynamicKmer{A}) where {A} = cmp(a.x, b.x)

@inline function Base.getindex(x::DynamicKmer{A}, idx::AbstractUnitRange{<:Integer}) where {A}
    @boundscheck checkbounds(x, idx)
    bps = BioSequences.bits_per_symbol(x)
    len = length(idx)
    u = left_shift(x.x, (first(idx) - 1) * bps)
    U = BioSequences.encoded_data_eltype(typeof(x))
    u &= top_mask(U, len * bps)
    return _new_dynamic_kmer(A, u | (len % U))
end

function BioSequences.complement(x::DynamicKmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    A = typeof(Alphabet(x))
    return _new_dynamic_kmer(A, x.x ⊻ coding_mask(x))
end

function BioSequences.complement(x::DynamicKmer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}})
    A = typeof(Alphabet(x))
    u = BioSequences.complement_bitpar(x.x, A())
    u &= coding_mask(x)
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

function Base.reverse(x::DynamicKmer{A}) where {A}
    Bps = BioSequences.BitsPerSymbol(A())
    u = BioSequences.reversebits(x.x, Bps)
    u = left_shift(u, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

function BioSequences.reverse_complement(x::DynamicKmer{<:NucleicAcidAlphabet})
    return reverse(complement(x))
end

BioSequences.iscanonical(x::DynamicKmer) = x <= reverse_complement(x)

# _n_gc

# Counting

# Construction utils!!!

# As integer, from integer. Should we guarantee the encoding scheme?
