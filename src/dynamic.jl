"""
    DynamicKmer{A <: Alphabet, U <: Unsigned} <: BioSequence{A}

Dynamic kmers are immutable, bitstype `BioSequence`s similar to `Kmer`s.
However, unlike the `Kmer` type, the length of a dynamic kmer is a run time
value, and not a compile time value.

Dynamic kmers are slightly less efficient than regular kmers.
They are useful when a workload includes kmers of varying sizes, where the
length specialization of the `Kmer` type would cause excessive compilation
and/or type instability.

See also: [`DynamicDNAKmer`](@ref), [`Kmer`](@ref)

# Examples
```jldoctest
julia> m = DynamicRNAKmer{UInt32}(rna"AUGCUGA")
7nt RNA Sequence:
AUGCUGA

julia> reverse_complement(m)
7nt RNA Sequence:
UCAGCAU
```
"""
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

utype(::Type{<:DynamicKmer{A, U}}) where {A, U} = U

"Alias for DynamicKmer{DNAAlphabet{2}, <:Unsigned}"
const DynamicDNAKmer{U} = (DynamicKmer{DNAAlphabet{2}, U} where {U <: Unsigned})

"Alias for DynamicKmer{RNAAlphabet{2}, <:Unsigned}"
const DynamicRNAKmer{U} = (DynamicKmer{RNAAlphabet{2}, U} where {U <: Unsigned})

"Alias for DynamicKmer{AminoAcidAlphabet, <:Unsigned}"
const DynamicAAKmer{U} = (DynamicKmer{AminoAcidAlphabet, U} where {U <: Unsigned})

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
    return iszero(x.x) ? x.x : top_mask(typeof(x), coding_bits(x))
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
    check_kmer(Kmer{A, K, N})
    length(x) == K || error("Must construct kmer from length K DynamicKmer")
    # This is now a compile time constant
    noncoding = 8 * sizeof(x) - BioSequences.bits_per_symbol(x) * K
    return if N == 0
        Kmer{A, K, N}(unsafe, ())
    elseif N == 1
        u = right_shift(x.x, noncoding) % UInt
        Kmer{A, K, N}(unsafe, (u,))
    else
        u = right_shift(x.x, noncoding)
        Nu = div(sizeof(x), sizeof(UInt))
        B = sizeof(x) * 8
        t = ntuple(Nu) do i
            right_shift(u, B - i * 8 * sizeof(UInt)) % UInt
        end
        Kmer{A, K, N}(unsafe, t)
    end
end

const HASH_MASK = 0x6ff6e9f0462d5162 % UInt

Base.copy(x::DynamicKmer) = x
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

function BioSequences._n_gc(x::DynamicKmer{<:NucleicAcidAlphabet})
    u = x.x & ~length_mask(typeof(x))
    return BioSequences.gc_bitcount(u, Alphabet(x))
end

function Kmers.as_integer(x::DynamicKmer)
    shift = (8 * sizeof(x) - coding_bits(x))
    return right_shift(x.x, shift)
end

function Kmers.from_integer(
        T::Type{DynamicKmer{A, U}}, x::U, len::Int
    ) where {A <: Alphabet, U <: Unsigned}
    if (len % UInt) > (capacity(T) % UInt)
        error("Length too large for dynamic kmer")
    end
    bps = BioSequences.bits_per_symbol(A())
    shift = 8 * sizeof(U) - len * bps
    u = left_shift(x, shift)
    return _new_dynamic_kmer(A, u | (len % U))
end

## More construction utils
function DynamicKmer{T1, U}(x::DynamicKmer{T2, U}) where {
        B,
        T1 <: NucleicAcidAlphabet{B},
        T2 <: NucleicAcidAlphabet{B},
        U <: Unsigned,
    }
    return _new_dynamic_kmer(T1, x.x)
end

function DynamicKmer{T1}(x::DynamicKmer{T2}) where {
        B,
        T1 <: NucleicAcidAlphabet{B},
        T2 <: NucleicAcidAlphabet{B},
    }
    return _new_dynamic_kmer(T1, x.x)
end

# Constructor dispatchtes to RecodingScheme
function DynamicKmer{A, U}(x) where {A <: Alphabet, U <: Unsigned}
    return build_dynamic_kmer(RecodingScheme(A(), typeof(x)), DynamicKmer{A, U}, x)
end

# Generic fallback for arbitrary iterables
function build_dynamic_kmer(::RecodingScheme, ::Type{T}, x) where {T}
    bps = BioSequences.bits_per_symbol(T)
    shift = 8 * sizeof(T)
    U = utype(T)
    u = zero(U)
    cap = capacity(T)
    len = 0
    A = Alphabet(T)
    for i in x
        len += 1
        shift -= bps
        len > cap && error("Iterator size exceeds maximum capacity of dynamic kmer")
        enc = BioSequences.encode(A, convert(eltype(T), i)) % U
        u |= left_shift(enc, shift)
    end
    return _new_dynamic_kmer(typeof(A), (len % U) | u)
end

# Here, we can extract the encoding directly
function build_dynamic_kmer(::Copyable, ::Type{T}, x::BioSequence) where {T}
    len = length(x)
    len > capacity(T) && error("Iterator size exceeds maximum capacity of dynamic kmer")
    bps = BioSequences.bits_per_symbol(T)
    shift = 8 * sizeof(T)
    U = utype(T)
    u = zero(U)
    A = Alphabet(T)
    for i in eachindex(x)
        shift -= bps
        enc = BioSequences.extract_encoded_element(x, i) % U
        u |= left_shift(enc, shift)
    end
    return _new_dynamic_kmer(typeof(A), (len % U) | u)
end

@inline function build_dynamic_kmer(
        R::AsciiEncode,
        ::Type{T},
        s::Union{String, SubString{String}},
    ) where {T}
    return build_dynamic_kmer(R, T, codeunits(s))
end

@inline function build_dynamic_kmer(::AsciiEncode, ::Type{T}, x::AbstractVector{UInt8}) where {T}
    len = length(x)
    len > capacity(T) && error("Iterator size exceeds maximum capacity of dynamic kmer")
    U = utype(T)
    u = zero(U)
    shift = 8 * sizeof(T)
    bps = BioSequences.bits_per_symbol(T)
    A = Alphabet(T)
    for i in eachindex(x)
        byte = x[i]
        encoding = BioSequences.ascii_encode(A, byte)
        if encoding > 0x7f
            throw(BioSequences.EncodeError(A, byte))
        end
        shift -= bps
        u |= left_shift(encoding % U, shift)
    end
    return _new_dynamic_kmer(typeof(A), (len % U) | u)
end

@inline function build_dynamic_kmer(::TwoToFour, ::Type{T}, x::BioSequence) where {T}
    len = length(x)
    len > capacity(T) && error("Iterator size exceeds maximum capacity of dynamic kmer")
    U = utype(T)
    u = zero(U)
    shift = 8 * sizeof(T)
    for i in eachindex(x)
        shift -= 4
        encoding = left_shift(one(U), BioSequences.extract_encoded_element(x, i) % U)
        u |= left_shift(encoding % U, shift)
    end
    return _new_dynamic_kmer(typeof(Alphabet(T)), (len % U) | u)
end

@inline function build_dynamic_kmer(::FourToTwo, ::Type{T}, x::BioSequence) where {T}
    len = length(x)
    len > capacity(T) && error("Iterator size exceeds maximum capacity of dynamic kmer")
    U = utype(T)
    u = zero(U)
    shift = 8 * sizeof(T)
    A = Alphabet(T)
    for i in eachindex(x)
        shift -= 2
        encoding = BioSequences.extract_encoded_element(x, i) % U
        isone(count_ones(encoding)) || throw_uncertain(A, eltype(x), encoding)
        u |= left_shift(trailing_zeros(encoding) % U, shift)
    end
    return _new_dynamic_kmer(typeof(A), (len % U) | u)
end

"""
    shift_encoding(x::DynamicKmer{A, U}, encoding::U)::typeof(x)

Add `encoding`, a valid encoding in the alphabet of the `x`,
and of the same integer type as that used in `x`,
to the end of dynamic kmer `x` and discarding the first symbol in `x`.

It is the user's responsibility to ensure that `encoding` is valid.

# Examples
```jldoctest
julia> enc = UInt32(0x0a); # encoding of DNA_Y in 4-bit alphabets

julia> kmer = DynamicKmer{DNAAlphabet{4}, UInt32}("TAGA");

julia> Kmers.shift_encoding(kmer, enc)
4nt DNA Sequence:
AGAY
```
"""
function shift_encoding(x::DynamicKmer{A, U}, encoding::U) where {A <: Alphabet, U <: Unsigned}
    mask = length_mask(typeof(x))
    u = x.x & ~mask
    u = left_shift(u, BioSequences.bits_per_symbol(x))
    u |= left_shift(encoding, noncoding_bits(x))
    u | (x.x & mask)
    _new_dynamic_kmer(A, u | (x.x & mask))
end