"""
    Oligo{A <: Alphabet, U <: Unsigned} <: BioSequence{A}

Dynamic kmers are immutable, bitstype `BioSequence`s similar to `Kmer`s.
However, unlike the `Kmer` type, the length of a dynamic kmer is a run time
value, and not a compile time value.

Dynamic kmers types have a maximum number of symbols they can store,
see [`capacity`](@ref) for details.

Dynamic kmers are slightly less efficient than regular kmers.
They are useful when a workload includes kmers of varying sizes, where the
length specialization of the `Kmer` type would cause excessive compilation
and/or type instability.

See also: [`DNAOligo`](@ref), [`Kmer`](@ref)

# Examples
```jldoctest
julia> m = RNAOligo{UInt32}(rna"AUGCUGA")
7nt RNAOligo{UInt32}:
AUGCUGA

julia> reverse_complement(m)
7nt RNAOligo{UInt32}:
UCAGCAU

julia> DNAKmer{7}(m)
DNA 7-mer:
ATGCTGA
```
"""
struct Oligo{A <: Alphabet, U <: Unsigned} <: BioSequence{A}
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

Base.summary(x::Oligo{<:Union{DNAAlphabet, RNAAlphabet}}) = string(length(x), "nt ", typeof(x))
Base.summary(x::Oligo{AminoAcidAlphabet}) = string(length(x), "aa ", typeof(x))
Base.summary(x::Oligo) = string(length(x), "-symbol ", typeof(x))

function Base.show(io::IO, ::MIME"text/plain", s::Oligo)
    println(io, summary(s), ':')
    return print(io, s)
end

Base.empty(::Type{<:Oligo{A, U}}) where {A, U} = _new_dynamic_kmer(A, zero(U))

utype(::Type{<:Oligo{A, U}}) where {A, U} = U

"""
    DNAOligo{U}
    DNAOligo(source)

Alias for `Oligo{DNAAlphabet{2}, U}`.
When called without an explicit backing type, constructs a `DNAOligo{UInt64}`.
"""
const DNAOligo{U} = (Oligo{DNAAlphabet{2}, U} where {U <: Unsigned})

"""
    RNAOligo{U}
    RNAOligo(source)

Alias for `Oligo{RNAAlphabet{2}, U}`.
When called without an explicit backing type, constructs an `RNAOligo{UInt64}`.
"""
const RNAOligo{U} = (Oligo{RNAAlphabet{2}, U} where {U <: Unsigned})

"""
    AAOligo{U}
    AAOligo(source)

Alias for `Oligo{AminoAcidAlphabet, U}`.
When called without an explicit backing type, constructs an `AAOligo{UInt64}`.
"""
const AAOligo{U} = (Oligo{AminoAcidAlphabet, U} where {U <: Unsigned})

Base.@constprop :aggressive Base.@assume_effects :foldable function max_coding_bits(
        ::Type{Oligo{A, U}}
    ) where {A, U}
    bps = BioSequences.bits_per_symbol(A())
    iszero(bps) && return 0
    n_bits = 8 * sizeof(U)
    for capacity in div(n_bits, bps):-1:0
        len_bits = 8 * sizeof(capacity) - leading_zeros(capacity)
        coding_bits = bps * capacity
        len_bits + coding_bits ≤ n_bits && return coding_bits
    end
    0
end

@inline function length_mask(T::Type{<:Oligo})
    U = BioSequences.encoded_data_eltype(T)
    return one(U) << (8 * sizeof(U) - max_coding_bits(T)) - one(U)
end

@inline function top_mask(::Type{U}, len::Integer) where {U <: Unsigned}
    return left_shift(typemax(U), (8 * sizeof(U) - len))
end

@inline function top_mask(T::Type{<:Oligo}, len::Integer)
    return top_mask(BioSequences.encoded_data_eltype(T), len)
end

@inline function coding_bits(x::Oligo)
    return BioSequences.bits_per_symbol(x) * length(x)
end

@inline function noncoding_bits(x::Oligo)
    return 8 * sizeof(x) - coding_bits(x)
end

@inline function coding_mask(x::Oligo)
    return iszero(x.x) ? x.x : top_mask(typeof(x), coding_bits(x))
end

"""
    capacity(T::Type{<:Oligo{A, U}})::Int
    capacity(x::Oligo)::Int

Compute the maximum number of symbols that an instance of the concrete
type `T` can contain.
This computation is a compile time constant and so should not take
any runtime computation.

The value of this number is not guaranteed to be stable across versions of Kmers.

If `B` is the number of bits per symbol of `A`, the answer is `clamp(typemax(U), Int)`,
else the answer is a number in `0:div(8 * sizeof(U), B)`

# Examples
```jldoctest
julia> capacity(DNAOligo{UInt32})
14

julia> capacity(AAOligo{UInt8})
0

julia> capacity(RNAOligo) # NB: UnionAll type
ERROR: MethodError: no method matching capacity(::Type{RNAOligo})
[...]
```
"""
Base.@constprop :aggressive Base.@assume_effects :foldable function capacity(
        T::Type{<:Oligo{A, U}}
    ) where {A, U}
    bps = BioSequences.bits_per_symbol(A())
    return if iszero(bps)
        # Clamp has a bug in Julia 1.10 and below which causes a normal `clamp`
        # call to overflow if typemax(U) > typemax(Int).
        # Therefore, use a manual implementation.
        if typemax(U) > typemax(Int)
            typemax(Int)
        else
            Int(typemax(U))::Int
        end
    else
        div(max_coding_bits(T), bps)
    end
end

@inline capacity(x::Oligo) = capacity(typeof(x))

BioSequences.encoded_data_eltype(::Type{Oligo{A, U}}) where {A, U} = U

function BioSequences.extract_encoded_element(x::Oligo, i::Integer)
    bps = BioSequences.bits_per_symbol(x)
    shift = 8 * sizeof(x) - (i * bps)
    u = right_shift(x.x, shift)
    mask = one(x.x) << (bps) - one(x.x)
    return u & mask
end

Base.length(x::Oligo) = (x.x & length_mask(typeof(x))) % Int
Base.isempty(x::Oligo) = iszero(x.x)

function Kmer{A, K}(x::Oligo{A}) where {A <: Alphabet, K}
    return @inline derive_type(Kmer{A, K})(x)
end

# This kmer construction and the one below is efficient, as kmers and dynamic kmers
# share a very similar encoding scheme.
function Kmer{A, K, N}(x::Oligo{A}) where {A <: Alphabet, K, N}
    check_kmer(Kmer{A, K, N})
    length(x) == K || error("Must construct kmer from length K Oligo")
    return from_integer(Kmer{A, K, N}, as_integer(x))
end

function Kmer{A1, K, N}(
        x::Oligo{A2}
    ) where {M, A1 <: NucleicAcidAlphabet{M}, A2 <: NucleicAcidAlphabet{M}, K, N}
    check_kmer(Kmer{A1, K, N})
    length(x) == K || error("Must construct kmer from length K Oligo")
    return from_integer(Kmer{A1, K, N}, as_integer(x))
end

const HASH_MASK = 0x6ff6e9f0462d5162 % UInt

Base.copy(x::Oligo) = x
Base.hash(x::Oligo, h::UInt) = hash(x.x, h ⊻ HASH_MASK)
fx_hash(x::Oligo, h::UInt) = ((bitrotate(h, 5) ⊻ x.x) % UInt) * FX_CONSTANT

Base.:(==)(a::Oligo{A, U}, b::Oligo{A, U}) where {A <: Alphabet, U <: Unsigned} = a === b

function Base.:(==)(
        a::Oligo{<:NucleicAcidAlphabet{N}, U},
        b::Oligo{<:NucleicAcidAlphabet{N}, U},
    ) where {N, U <: Unsigned}
    return a.x == b.x
end

# Oligos are compared only when their backing representations are compatible.
# Comparing another Oligo or BioSequence would otherwise fall back to decoded
# BioSymbol comparison, which cannot share this type's representation-based hash.
Base.:(==)(a::Oligo, b::Oligo) = throw(MethodError(==, (a, b)))
Base.:(==)(a::Oligo, b::BioSequence) = throw(MethodError(==, (a, b)))
Base.:(==)(a::BioSequence, b::Oligo) = throw(MethodError(==, (a, b)))

# These resolve the intersection with Kmer's corresponding BioSequence fallbacks.
Base.:(==)(a::Oligo, b::Kmer) = throw(MethodError(==, (a, b)))
Base.:(==)(a::Kmer, b::Oligo) = throw(MethodError(==, (a, b)))

function Base.isless(
        a::Oligo{A, U}, b::Oligo{A, U}
    ) where {A <: Alphabet, U <: Unsigned}
    return isless(a.x, b.x)
end

function Base.cmp(
        a::Oligo{A, U}, b::Oligo{A, U}
    ) where {A <: Alphabet, U <: Unsigned}
    return cmp(a.x, b.x)
end

function Base.isless(
        a::Oligo{<:NucleicAcidAlphabet{N}, U},
        b::Oligo{<:NucleicAcidAlphabet{N}, U},
    ) where {N, U <: Unsigned}
    return isless(a.x, b.x)
end

function Base.cmp(
        a::Oligo{<:NucleicAcidAlphabet{N}, U},
        b::Oligo{<:NucleicAcidAlphabet{N}, U},
    ) where {N, U <: Unsigned}
    return cmp(a.x, b.x)
end

Base.isless(a::Oligo, b::Oligo) = throw(MethodError(isless, (a, b)))
Base.cmp(a::Oligo, b::Oligo) = throw(MethodError(cmp, (a, b)))

"""
    widen(mer::Oligo)
    widen(::Type{T}) where {T <: Oligo}

Given an `Oligo` type, return the `Oligo` type with the same alphabet
and with a larger backing integer type. E.g. for a `Oligo{A, UInt32}`,
return `Oligo{A, UInt64}`.
If there are no defined bit-integers wider than the one backing the input type,
throw a `MethodError`.

For the input integer types `UInt8`, `UInt16`, `UInt32`, `UInt64` and `UInt128`,
the widened integer type will be next in that sequence (erroring for `UInt128`
as input).
When the `BitIntegers` package is loaded, this sequence extends to
`UInt256`, `UInt512` and `UInt1024`.

When given an instance `x::Oligo`, construct `widen(typeof(x))(x)`.

# Examples
```
julia> m1 = DNAOligo{UInt16}("TAGCTAG"); m2 = widen(m1);

julia> m2 === DNAOligo{UInt32}(m1)
true
```
"""
Base.widen(x::Oligo) = widen(typeof(x))(x)

function Base.widen(::Type{Oligo{A, U}}) where {A, U}
    return Oligo{A, widen_bitint(U)}
end

# We use this instead of Base.widen, because widen(UInt128) == BigInt,
# and we don't back Oligo by BigInt.
widen_bitint(T::Union{Type{UInt8}, Type{UInt16}, Type{UInt32}, Type{UInt64}}) = widen(T)

@inline function Base.getindex(x::Oligo{A}, idx::AbstractUnitRange{<:Integer}) where {A}
    @boundscheck if isempty(idx)
        i = first(idx)
        if i < firstindex(x) || (i > lastindex(x) && i - one(i) != lastindex(x))
            boundserror(x, idx)
        end
    else
        checkbounds(x, first(idx))
        checkbounds(x, last(idx))
    end
    isempty(idx) && return empty(typeof(x))
    bps = BioSequences.bits_per_symbol(x)
    len = length(idx)
    u = left_shift(x.x, (first(idx) - 1) * bps)
    U = BioSequences.encoded_data_eltype(typeof(x))
    u &= top_mask(U, len * bps)
    return _new_dynamic_kmer(A, u | (len % U))
end

@inline function Base.getindex(
        x::Oligo{A, U}, indices::AbstractVector{Bool}
    ) where {A, U <: Unsigned}
    @boundscheck checkbounds(x, indices)
    bps = BioSequences.bits_per_symbol(x)
    source = x.x
    u = zero(U)
    len = 0
    for selected in indices
        if selected
            len += 1
            encoding = if iszero(bps)
                zero(U)
            else
                right_shift(source, 8 * sizeof(U) - bps)
            end
            u = left_shift(u, bps) | encoding
        end
        source = left_shift(source, bps)
    end
    u = left_shift(u, 8 * sizeof(U) - len * bps)
    return _new_dynamic_kmer(A, u | (len % U))
end

function Base.getindex(
        x::Oligo{A, U}, indices::AbstractVector{<:Integer}
    ) where {A, U <: Unsigned}
    len = length(indices)
    if len > capacity(typeof(x))
        boundserror(x, indices)
    end
    bps = BioSequences.bits_per_symbol(x)
    shift = 8 * sizeof(U)
    u = zero(U)
    for i in indices
        @boundscheck checkbounds(x, i)
        shift -= bps
        encoding = BioSequences.extract_encoded_element(x, i) % U
        u |= left_shift(encoding, shift)
    end
    return _new_dynamic_kmer(A, u | (len % U))
end

function BioSequences.complement(x::Oligo{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    A = typeof(Alphabet(x))
    return _new_dynamic_kmer(A, x.x ⊻ coding_mask(x))
end

function BioSequences.complement(x::Oligo{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}})
    A = typeof(Alphabet(x))
    u = BioSequences.complement_bitpar(x.x, A())
    u &= coding_mask(x)
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

# Generic fallback for nonstandard nucleic-acid alphabets. The specialized
# two-bit and four-bit methods above avoid decoding and re-encoding.
function BioSequences.complement(x::Oligo{<:NucleicAcidAlphabet})
    return typeof(x)((complement(i) for i in x))
end

function Base.reverse(x::Oligo{A}) where {A}
    Bps = BioSequences.BitsPerSymbol(A())
    u = BioSequences.reversebits(x.x, Bps)
    u = left_shift(u, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

function BioSequences.reverse_complement(x::Oligo{<:NucleicAcidAlphabet})
    return reverse(complement(x))
end

BioSequences.iscanonical(x::Oligo) = x <= reverse_complement(x)

# This is more efficient than the fallback because RC'ing is cheap
function BioSequences.canonical(x::Oligo{<:NucleicAcidAlphabet})
    rc = reverse_complement(x)
    return x < rc ? x : rc
end

function BioSequences._n_gc(x::Oligo{<:NucleicAcidAlphabet})
    u = x.x & ~length_mask(typeof(x))
    return BioSequences.gc_bitcount(u, Alphabet(x))
end

"""
    as_integer(x::Oligo{A, U})::U

Similar to `as_integer` for kmers, but is guaranteed to return a value of `U`,
and the number of coding bits is known at runtime.
"""
function Kmers.as_integer(x::Oligo)
    shift = (8 * sizeof(x) - coding_bits(x))
    return right_shift(x.x, shift)
end

"""
    from_integer(T::Type{<:Oligo{A, U}}, u::Unsigned, len::Int)::T

Similar to `from_integer` for `Kmer`, but the length of the resulting `Oligo`
must be passed as an argument. Will error if `len` is larger than the maximal size
supported by `T`.

The input may be any unsigned integer type. Only the lowest coding bits are used,
so round trips are independent of the width of the input integer when it can hold
all coding bits.

If `u` is obtained from a `Oligo` with a length different from `len`,
the resulting `Oligo` is reproducible, but not correct and may change between
versions.

# Examples
```jldoctest
julia> d = DNAOligo{UInt32}(dna"TAGTGCTGTAGGC")
13nt DNAOligo{UInt32}:
TAGTGCTGTAGGC

julia> u = as_integer(d);

julia> from_integer(typeof(d), u, 13) === d
true

julia> from_integer(typeof(d), u, 12) == d
false
```
"""
function Kmers.from_integer(
        T::Type{Oligo{A, U}}, x::Unsigned, len::Int
    ) where {A <: Alphabet, U <: Unsigned}
    if (len % UInt) > (capacity(T) % UInt)
        error("Length too large for dynamic kmer")
    end
    bps = BioSequences.bits_per_symbol(A())
    iszero(bps) && return _new_dynamic_kmer(A, len % U)
    iszero(len) && return _new_dynamic_kmer(A, zero(U))
    shift = 8 * sizeof(U) - len * bps
    u = left_shift(x % U, shift)
    return _new_dynamic_kmer(A, u | (len % U))
end

################################################
# Unsafe extract
################################################

@inline function extract_elements(
        recoding::RecodingScheme,
        ::Type{T},
        source,
        from::Int,
        len::Int,
    ) where {A, U, T <: Oligo{A, U}}
    u = zero(U)
    alphabet = A()
    bps = BioSequences.bits_per_symbol(alphabet)
    for i in 0:(len - 1)
        encoding = recode_element(recoding, alphabet, source, from + i, U)
        u = left_shift(u, bps) | encoding
    end
    u = left_shift(u, 8 * sizeof(U) - bps * len)
    return _new_dynamic_kmer(A, u | (len % U))
end

@inline function unsafe_extract(
        recoding::RecodingScheme,
        ::Type{T},
        source,
        from::Int,
        len::Int,
    ) where {A, U, T <: Oligo{A, U}}
    return extract_elements(recoding, T, source, from, len)
end

@inline function packed_bits(
        ::Copyable,
        ::Type{U},
        ::Alphabet,
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
    ) where {U <: Unsigned}
    bps = BioSequences.BitsPerSymbol(source)
    if sizeof(U) <= sizeof(UInt64)
        return load_packed_bits(UInt64, source, from, len, bps) % U
    end
    return load_packed_bits(U, source, from, len, bps)
end

@inline function packed_bits(
        ::TwoToFour,
        ::Type{U},
        ::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{2}},
        from::Int,
        len::Int,
    ) where {U <: Unsigned}
    bits = zero(U)
    offset = 0
    while offset < len
        chunk_len = min(16, len - offset)
        chunk = load_packed_bits(
            UInt64,
            source,
            from + offset,
            chunk_len,
            BioSequences.BitsPerSymbol{2}(),
        )
        encoding = BioSequences.two_to_four_bits(chunk % UInt32) &
            BioSequences.bitmask(UInt64, 4 * chunk_len)
        bits |= left_shift(encoding % U, 4 * offset)
        offset += chunk_len
    end
    return bits
end

@inline function packed_bits(
        ::FourToTwo,
        ::Type{U},
        alphabet::Alphabet,
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        from::Int,
        len::Int,
    ) where {U <: Unsigned}
    bits = zero(U)
    offset = 0
    while offset < len
        chunk_len = min(16, len - offset)
        encoding = four_to_two_half(alphabet, source, from + offset, chunk_len)
        bits |= left_shift(encoding % U, 2 * offset)
        offset += chunk_len
    end
    return bits
end

@inline function reverse_packed_bits(bits::Unsigned, bps::BioSequences.BitsPerSymbol)
    return BioSequences.reversebits(bits, bps)
end

# BioSequences delegates one-bit reversal to `bitreverse`, which does not support
# wide BitIntegers. Reverse those one UInt64 at a time instead.
@inline function reverse_packed_bits(
        bits::U,
        ::BioSequences.BitsPerSymbol{1},
    ) where {U <: Unsigned}
    if sizeof(U) <= sizeof(UInt64)
        shift = 8 * (sizeof(UInt64) - sizeof(U))
        return right_shift(Base.bitreverse(bits % UInt64), shift) % U
    end

    width = 8 * sizeof(U)
    reversed = zero(U)
    for offset in 0:64:(width - 1)
        nbits = min(64, width - offset)
        chunk = right_shift(bits, offset) % UInt64
        chunk = right_shift(Base.bitreverse(chunk), 64 - nbits)
        reversed |= left_shift(chunk % U, width - offset - nbits)
    end
    return reversed
end

@inline function extract_packed(
        recoding::Union{Copyable, TwoToFour, FourToTwo},
        ::Type{T},
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
        bps::BioSequences.BitsPerSymbol{B},
    ) where {A, U, T <: Oligo{A, U}, B}
    iszero(B) && return _new_dynamic_kmer(A, len % U)
    iszero(len) && return _new_dynamic_kmer(A, zero(U))

    alphabet = A()
    bits = packed_bits(recoding, U, alphabet, source, from, len)
    u = reverse_packed_bits(bits, bps)
    return _new_dynamic_kmer(A, u | (len % U))
end

@inline function unsafe_extract(
        recoding::Copyable,
        ::Type{T},
        source::BioSequences.SeqOrView,
        from::Int,
        len::Int,
    ) where {A, U, T <: Oligo{A, U}}
    bps = BioSequences.BitsPerSymbol(A())
    if bps == BioSequences.BitsPerSymbol(source) && !iszero(BioSequences.bits_per_symbol(bps))
        return extract_packed(recoding, T, source, from, len, bps)
    end
    return extract_elements(recoding, T, source, from, len)
end

@inline function unsafe_extract(
        recoding::TwoToFour,
        ::Type{T},
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{2}},
        from::Int,
        len::Int,
    ) where {A <: NucleicAcidAlphabet{4}, U, T <: Oligo{A, U}}
    return extract_packed(recoding, T, source, from, len, BioSequences.BitsPerSymbol(A()))
end

@inline function unsafe_extract(
        recoding::FourToTwo,
        ::Type{T},
        source::BioSequences.SeqOrView{<:NucleicAcidAlphabet{4}},
        from::Int,
        len::Int,
    ) where {A <: NucleicAcidAlphabet{2}, U, T <: Oligo{A, U}}
    return extract_packed(recoding, T, source, from, len, BioSequences.BitsPerSymbol(A()))
end

## More construction utils
function Oligo{T1, U}(x::Oligo{T2, U}) where {
        B,
        T1 <: NucleicAcidAlphabet{B},
        T2 <: NucleicAcidAlphabet{B},
        U <: Unsigned,
    }
    return _new_dynamic_kmer(T1, x.x)
end

"""
    Oligo{A}(source)
    DNAOligo(source)
    RNAOligo(source)
    AAOligo(source)

Construct an `Oligo` with alphabet `A` and a `UInt64` backing integer from `source`.
To choose another backing type, use `Oligo{A, U}(source)` or the
corresponding `*Oligo{U}` alias.
"""
@inline function Oligo{A}(x) where {A <: Alphabet}
    return Oligo{A, UInt64}(x)
end

# Constructor dispatches to RecodingScheme
@inline function Oligo{A, U}(x) where {A <: Alphabet, U <: Unsigned}
    return build_dynamic_kmer(RecodingScheme(A(), typeof(x)), Oligo{A, U}, x)
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

# Indexed BioSequences can all use the public runtime-length extraction API.
@inline function build_dynamic_kmer(R::RecodingScheme, ::Type{T}, x::BioSequence) where {T}
    len = length(x)
    len > capacity(T) && error("Iterator size exceeds maximum capacity of dynamic kmer")
    return unsafe_extract(R, T, x, firstindex(x), len)
end

# More efficient, since the internal representation is even closer.
# This function is huge, but a lot of it is compile time work, so force inline it.
@inline function build_dynamic_kmer(::Copyable, ::Type{T}, x::Kmer) where {T}
    len = length(x)
    len > capacity(T) && error("Kmer size exceeds maximum capacity of dynamic kmer")
    A = Alphabet(T)
    U = utype(T)
    bps = BioSequences.bits_per_symbol(A)
    # If no BPS, binary representation of dynamic kmer is no coding bits,
    # and then simply the length in the lower bits.
    iszero(bps) && return _new_dynamic_kmer(typeof(A), len % U)
    isempty(x) &&return _new_dynamic_kmer(typeof(A), zero(U))
    tup = BioSequences.encoded_data(x)
    # Tuple length has to be at least one, since otherwise kmer would be empty,
    # and branch above would have been taken

    # Kmers store the partially filled head word in its lower bits, whereas an
    # Oligo stores all coding bits at the top of its backing integer. Align
    # the head word relative to the target width, rather than the source UInt
    # width, so this also works when U is wider than UInt.
    u = zero(U)
    shift = 8 * sizeof(U) - (64 - bits_unused(typeof(x)))
    for word in tup
        u |= left_shift(word % U, shift)
        shift -= 64
    end

    # Add in length
    u |= len % U
    return _new_dynamic_kmer(typeof(A), u)
end

function build_dynamic_kmer(::Copyable, ::Type{T}, x::Oligo) where {T}
    u = @inline switch_backing_encoding(utype(T), x)
    A = Alphabet(T)
    return _new_dynamic_kmer(typeof(A), u)
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
    return unsafe_extract(AsciiEncode(), T, x, firstindex(x), len)
end

# Switch encoding data of `x` to `T`. Error if it doesn't fit.
function switch_backing_encoding(T::Type{<:Unsigned}, x::Oligo{A, U}) where {A, U}
    T == U && return x
    return if sizeof(T) < sizeof(x)
        narrow_to(T, x)
    else
        widen_to(T, x)
    end
end

# Create a Oligo{A, T} containig the same sequence as `x`, efficiently,
# or error if `x` does not fit in that type.
function narrow_to(T::Type{<:Unsigned}, x::Oligo{A, U}) where {A, U}
    newT = Oligo{A, T}
    if length(x) > capacity(newT)
        error("Dynamic Kmer do not fit into integer size")
    end
    # Remove length from encoding
    mask = length_mask(typeof(x))
    u = x.x & ~mask

    # Shift down to new location in smaller integer
    shift = 8 * (sizeof(U) - sizeof(T))
    u = right_shift(u, shift)

    # Add length back and return
    u |= (x.x & mask)
    return u % T
end

# Create a Oligo{A, T} containig the same sequence as `x`, efficiently
function widen_to(T::Type{<:Unsigned}, x::Oligo{A, U}) where {A, U}
    # Remove length from encoding
    mask = length_mask(typeof(x))
    u = (x.x & ~mask) % T

    # Shift up to new location in larger integer
    shift = 8 * (sizeof(T) - sizeof(U))
    u = left_shift(u, shift)

    # Add length back and return
    u |= (x.x & mask)
    return u % T
end

"""
    shift_encoding(x::Oligo{A, U}, encoding::U)::typeof(x)

Add `encoding`, a valid encoding in the alphabet of the `x`,
and of the same integer type as that used in `x`,
to the end of dynamic kmer `x` and discarding the first symbol in `x`.

It is the user's responsibility to ensure that `encoding` is valid.

# Examples
```jldoctest
julia> enc = UInt32(0x0a); # encoding of DNA_Y in 4-bit alphabets

julia> kmer = Oligo{DNAAlphabet{4}, UInt32}("TAGA");

julia> Kmers.shift_encoding(kmer, enc)
4nt Oligo{DNAAlphabet{4}, UInt32}:
AGAY
```
"""
@inline function shift_encoding(
        x::Oligo{A, U}, encoding::U
    ) where {A <: Alphabet, U <: Unsigned}
    isempty(x) && return x
    mask = length_mask(typeof(x))
    u = x.x & ~mask
    bps = BioSequences.bits_per_symbol(x)
    iszero(bps) && return x
    u = left_shift(u, bps)
    u |= left_shift(encoding, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & mask))
end

"""
    shift(x::T, s)::T where {T <: Oligo}

Push `s` onto the end of `x` and discard the first symbol. The alphabet,
backing integer type, and length are preserved. The argument `s` is converted
to the element type of `x` before it is encoded.

Shifting an empty `Oligo` returns it unchanged.

See also: [`shift_first`](@ref), [`push`](@ref)

# Examples
```jldoctest
julia> shift(dmer"TACC"d, RNA_A)
4nt DNAOligo{UInt64}:
ACCA

julia> shift(dmer"KWOP"a, 'L')
4aa AAOligo{UInt64}:
WOPL
```
"""
@inline function shift(x::Oligo{A, U}, s) where {A <: Alphabet, U <: Unsigned}
    E = eltype(x)
    symbol = convert(E, s)::E
    encoding = U(BioSequences.encode(A(), symbol))::U
    return shift_encoding(x, encoding)
end

"""
    shift_first(x::T, s)::T where {T <: Oligo}

Push `s` onto the start of `x` and discard the last symbol. The alphabet,
backing integer type, and length are preserved. The argument `s` is converted
to the element type of `x` before it is encoded.

Shifting an empty `Oligo` returns it unchanged.

See also: [`shift`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> shift_first(dmer"TACC"d, RNA_A)
4nt DNAOligo{UInt64}:
ATAC

julia> shift_first(dmer"KWOP"a, 'L')
4aa AAOligo{UInt64}:
LKWO
```
"""
@inline function shift_first(x::Oligo{A, U}, s) where {A <: Alphabet, U <: Unsigned}
    E = eltype(x)
    symbol = convert(E, s)::E
    encoding = U(BioSequences.encode(A(), symbol))::U
    return shift_first_encoding(x, encoding)
end

"""
    shift_first_encoding(x::Oligo{A, U}, encoding::U)::typeof(x)

Add `encoding`, a valid encoding in the alphabet of `x`, to the start of `x`
and discard its last symbol. It is the caller's responsibility to ensure that
`encoding` is valid.
"""
@inline function shift_first_encoding(
        x::Oligo{A, U}, encoding::U
    ) where {A <: Alphabet, U <: Unsigned}
    isempty(x) && return x
    bps = BioSequences.bits_per_symbol(x)
    iszero(bps) && return x
    mask = length_mask(typeof(x))
    u = right_shift(x.x & ~mask, bps)
    u &= coding_mask(x)
    u |= left_shift(encoding, 8 * sizeof(U) - bps)
    return _new_dynamic_kmer(A, u | (x.x & mask))
end

Base.adjoint(x::Oligo) = x

"""
    @dmer_str -> Oligo

Construct a `Oligo{A, UInt64}` from the given string. The macro must be used with a flag
after the input string, e.g. `d` in `dmer"TAG"d` or `a` in `dmer"PCW"a`, signifying
the alphabet of the dynamic kmer.
The flags `d = DNAAlphabet{2}`, `r = RNAAlphabet{2}` and `a = AminoAcidAlphabet`
are recognized.

# Examples
```jldoctest
julia> dmer"UGCUA"r
5nt RNAOligo{UInt64}:
UGCUA

julia> dmer"YDLLKKR"a
7aa AAOligo{UInt64}:
YDLLKKR

julia> dmer"TATTAGCA"d
8nt DNAOligo{UInt64}:
TATTAGCA
```
"""
macro dmer_str(seq, flag)
    trimmed = BioSequences.remove_newlines(seq)
    # Unlike @dna_str, we default to 2-bit alphabets, because kmers
    # by convention are usually 2-bit only
    return if flag == "dna" || flag == "d"
        DNAOligo{UInt64}(trimmed)
    elseif flag == "rna" || flag == "r"
        RNAOligo{UInt64}(trimmed)
    elseif flag == "aa" || flag == "a"
        AAOligo{UInt64}(trimmed)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

"""
    push(x::T, s)::T where {T <: Oligo}

Create a new `Oligo` of type `T` by adding the symbol `s` to the end of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `BoundsError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push_first`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DNAOligo{UInt64}:
TGTGCTGA

julia> d2 = push(d, 'G') # converts from Char to DNA
9nt DNAOligo{UInt64}:
TGTGCTGAG

julia> d == d2 # does not mutate immutable d
false

julia> push(dmer"RRKRLVD"a, AA_W)
ERROR: BoundsError: attempt to access AAOligo{UInt64} at index [8]
[...]
```
"""
function push(x::Oligo{A, U}, s) where {A, U}
    T = typeof(x)

    len = length(x)
    @boundscheck if len >= capacity(T)
        boundserror(x, capacity(T) + 1)
    end
    new_len = len + 1

    # Update the packed length only after its logical value has been validated.
    u = x.x + one(U)

    E = eltype(x)
    sT = convert(E, s)::E
    enc = U(BioSequences.encode(A(), sT))::U
    bps = BioSequences.bits_per_symbol(A())

    shift = (8 * sizeof(U)) - (bps * new_len)
    u |= left_shift(enc, shift)
    return _new_dynamic_kmer(A, u)
end

@noinline boundserror(x, i) = throw(BoundsError(x, i))

"""
    push_first(x::T, s)::T where {T <: Oligo}

Create a new `Oligo` of type `T` by adding the symbol `s` to the start of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `BoundsError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DNAOligo{UInt64}:
TGTGCTGA

julia> d2 = push_first(d, 'G') # converts from Char to DNA
9nt DNAOligo{UInt64}:
GTGTGCTGA

julia> d == d2 # does not mutate immutable d
false

julia> push_first(dmer"RRKRLVD"a, AA_W)
ERROR: BoundsError: attempt to access AAOligo{UInt64} at index [8]
[...]
```
"""
function push_first(x::Oligo{A, U}, s) where {A, U}
    T = typeof(x)

    len = length(x)
    @boundscheck if len >= capacity(T)
        boundserror(x, capacity(T) + 1)
    end
    new_len = len + 1


    E = eltype(x)
    sT = convert(E, s)::E
    enc = U(BioSequences.encode(A(), sT))::U
    bps = BioSequences.bits_per_symbol(A())

    mask = length_mask(T)
    # Remove length. Since we shift, the length would be garbled
    u = x.x & ~mask
    # Shift down to make room for symbol at head
    u >>= bps
    # Add in symbol at head
    shift = 8 * sizeof(U) - bps
    u |= enc << shift

    # Add in new length
    u |= new_len % U
    return _new_dynamic_kmer(A, u)
end

"""
    pop(x::Oligo{A, U})::Oligo{A, U}

Returns a new dynamic kmer with the last symbol of the input `x` removed.
Throws an `BoundsError` if `x` is empty.

See also: [`pop_first`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"EDEAVY"a
6aa AAOligo{UInt64}:
EDEAVY

julia> d2 = pop(d)
5aa AAOligo{UInt64}:
EDEAV

julia> d == d2
false

julia> pop(dmer""a)
ERROR: BoundsError: attempt to access AAOligo{UInt64} at index [0]
[...]
```
"""
function pop(x::Oligo{A, U}) where {A, U}
    isempty(x) && boundserror(x, 0)

    # Decrement length
    u = x.x
    u -= one(u)

    # Remove the symbol
    bps = BioSequences.bits_per_symbol(A())
    mask = U(1) << bps - U(1)
    shift = 8 * sizeof(U) - length(x) * bps
    u &= ~left_shift(mask, shift)
    return _new_dynamic_kmer(A, u)
end

"""
    pop_first(x::Oligo{A, U})::Oligo{A, U}

Returns a new dynamic kmer with the first symbol of the input `x` removed.
Throws an `BoundsError` if `x` is empty.

See also: [`pop`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"UGCGUAGCUA"r
10nt RNAOligo{UInt64}:
UGCGUAGCUA

julia> d2 = pop_first(d)
9nt RNAOligo{UInt64}:
GCGUAGCUA

julia> d == d2
false

julia> pop_first(dmer""r)
ERROR: BoundsError: attempt to access RNAOligo{UInt64} at index [0]
[...]
```
"""
function pop_first(x::Oligo{A, U}) where {A, U}
    isempty(x) && boundserror(x, 0)

    # Remove length, since we need to shift it to pop first,
    # and shifting would move the length bits
    mask = length_mask(Oligo{A, U})
    u = x.x & ~mask

    # Remove the symbol by shifting
    bps = BioSequences.bits_per_symbol(A())
    u <<= bps

    # Add in length back minus one
    u |= (x.x & mask) - 0x01
    return _new_dynamic_kmer(A, u)
end

@noinline throw_argumenterror(s::String) = throw(ArgumentError(s))

function Base.setindex(kmer::Oligo{A, U}, v, i::Integer) where {A, U}
    i = Int(i)::Int
    @boundscheck checkbounds(kmer, i)

    # Convert to the element type of A, then encode to a U
    E = eltype(typeof(kmer))
    sT = convert(E, v)::E
    enc = U(BioSequences.encode(A(), sT))::U
    bps = BioSequences.bits_per_symbol(A())
    iszero(bps) && return kmer

    # Zero out the bits that code ofr the i'th symbol in `kmer`
    shift = 8 * sizeof(U) - i * bps
    mask = U(1) << bps - U(1)
    u = kmer.x
    u &= ~left_shift(mask, shift)

    # Now add in the encoding bits at the right location and return kmer
    u |= left_shift(enc, shift)
    return _new_dynamic_kmer(A, u)
end

"""
    translate(
        seq::Oligo{<:Union{DNAAlphabet, RNAAlphabet}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true,
        alternative_start::Bool = false,
    )::AAOligo

Translate a nucleotide `Oligo` to a `AAOligo`.
The type of the result is the smallest `AAOligo`, which is statically known to
have a capacity large enough to hold the result.

If the result does not fit in the largest known `AAOligo`, throw an exception.
You can increase the largest known `AAOligo` by loading the package
BitIntegers.jl.

The arguments other than `seq` are identical to the method with `LongSequence`.

# Examples
```jldoctest
julia> d = Oligo{DNAAlphabet{4}, UInt64}("TGGCCCGATTGA");

julia> translate(dmer"TGGCCCGATTGA"d)
4aa AAOligo{UInt128}:
WPD*

julia> translate(DNAOligo{UInt32}("TGGCCCGATTGA"); alternative_start=true)
4aa AAOligo{UInt64}:
MPD*
```
"""
function BioSequences.translate(
        seq::Oligo{<:Union{DNAAlphabet, RNAAlphabet}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true,
        alternative_start::Bool = false,
    )::AAOligo

    # Check length of sequence is divisible by 3
    (aalen, remainder) = divrem(length(seq) % UInt, 3 % UInt)
    iszero(remainder) || throw_argumenterror("Dynamic kmer length not divisible by 3")

    # Build the encoding. If `alternative_start`, we need to set the starting
    # AA to AA_M, no matter what. So, we simply skip the first codon.
    U = get_matching_aaseq_utype(typeof(seq))
    u = zero(U)

    # Begin by shifting to top bits, or, if alternative_start, top bits
    # except the top 8 bits which should have the AA_M
    shift = 8 * sizeof(U) - 8 - 8 * alternative_start
    for i in (1 + alternative_start):(aalen % Int)
        aa = inbounds_aa_from(seq, code, 3i - 2, allow_ambiguous_codons)
        encoding = reinterpret(UInt8, aa) % U
        u |= left_shift(encoding, shift)
        shift -= 8
    end

    u |= aalen % U

    # If alternative_start, we manually add the AA_M (encoding 0x0c)
    # to the data at the top bits
    if alternative_start && !iszero(aalen)
        u |= left_shift(0x0c % U, 8 * sizeof(U) - 8)
    end

    return _new_dynamic_kmer(AminoAcidAlphabet, u)
end

@inline function inbounds_aa_from(
        seq::Oligo{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}},
        code::BioSequences.GeneticCode,
        i::Int,
        _::Bool,
    )
    a = BioSequences.extract_encoded_element(seq, i)
    b = BioSequences.extract_encoded_element(seq, i + 1)
    c = BioSequences.extract_encoded_element(seq, i + 2)
    codon = (a << 4) | (b << 2) | c
    return @inbounds code[codon % UInt64]
end

@inline function inbounds_aa_from(
        seq::Oligo{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}},
        code::BioSequences.GeneticCode,
        i::Int,
        allow_ambiguous_codons::Bool,
    )
    a = @inbounds reinterpret(RNA, seq[i])
    b = @inbounds reinterpret(RNA, seq[i + 1])
    c = @inbounds reinterpret(RNA, seq[i + 2])
    return if isgap(a) | isgap(b) | isgap(c)
        error("Cannot translate nucleotide sequences with gaps.")
    elseif iscertain(a) & iscertain(b) & iscertain(c)
        code[BioSequences.unambiguous_codon(a, b, c)]
    else
        BioSequences.try_translate_ambiguous_codon(code, a, b, c, allow_ambiguous_codons)
    end
end

function get_large_bitsize end

@inline Base.@constprop :aggressive Base.@assume_effects :foldable function get_matching_aaseq_utype(
        T::Type{<:Oligo{<:NucleicAcidAlphabet}}
    )
    max_aa = div(capacity(T) % UInt, UInt(3)) % Int

    return if max_aa < 2
        UInt16
    elseif max_aa < 4
        UInt32
    elseif max_aa < 8
        UInt64
    elseif max_aa < 16
        UInt128
    else
        if hasmethod(get_large_bitsize, Tuple{Val})
            get_large_bitsize(Val{max_aa}())
        else
            error(
                "AA oligo does not fit in 128 bits. " *
                    "Load package BitIntegers to access Unsigned types larger than UInt128."
            )
        end
    end
end
