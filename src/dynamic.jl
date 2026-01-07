"""
    DynamicKmer{A <: Alphabet, U <: Unsigned} <: BioSequence{A}

Dynamic kmers are immutable, bitstype `BioSequence`s similar to `Kmer`s.
However, unlike the `Kmer` type, the length of a dynamic kmer is a run time
value, and not a compile time value.

Dynamic kmers types have a maximum number of symbols they can store,
see [`capacity`](@ref) for details.

Dynamic kmers are slightly less efficient than regular kmers.
They are useful when a workload includes kmers of varying sizes, where the
length specialization of the `Kmer` type would cause excessive compilation
and/or type instability.

See also: [`DynamicDNAKmer`](@ref), [`Kmer`](@ref)

# Examples
```jldoctest
julia> m = DynamicRNAKmer{UInt32}(rna"AUGCUGA")
7nt DynamicRNAKmer{UInt32}:
AUGCUGA

julia> reverse_complement(m)
7nt DynamicRNAKmer{UInt32}:
UCAGCAU

julia> DNAKmer{7}(m)
DNA 7-mer:
ATGCTGA
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

Base.summary(x::DynamicKmer{<:Union{DNAAlphabet, RNAAlphabet}}) = string(length(x), "nt ", typeof(x))
Base.summary(x::DynamicKmer{AminoAcidAlphabet}) = string(length(x), "aa ", typeof(x))
Base.summary(x::DynamicKmer) = string(length(x), "-symbol ", typeof(x))

function Base.show(io::IO, ::MIME"text/plain", s::DynamicKmer)
    println(io, summary(s), ':')
    return print(io, s)
end

Base.empty(::Type{<:DynamicKmer{A, U}}) where {A, U} = _new_dynamic_kmer(A, zero(U))

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
    n_bits = 8 * sizeof(U)
    for capacity in div(n_bits, bps):-1:0
        len_bits = 8 * sizeof(capacity) - leading_zeros(capacity)
        coding_bits = bps * capacity
        len_bits + coding_bits ≤ n_bits && return coding_bits
    end
    0
end

@inline function length_mask(T::Type{<:DynamicKmer})
    U = BioSequences.encoded_data_eltype(T)
    return one(U) << (8 * sizeof(U) - max_coding_bits(T)) - one(U)
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

"""
    capacity(T::Type{<:DynamicKmer{A, U}})::Int

Compute the maximum number of symbols that an instance of the concrete
type `T` can contain.
This computation is a compile time constant and so should not take
any runtime computation.

The value of this number is not guaranteed to be stable across versions of Kmers.

If `B` is the number of bits per symbol of `A`, the answer is `clamp(typemax(U), Int)`,
else the answer is a number in `0:div(8 * sizeof(U), B)`

# Examples
```jldoctest
julia> capacity(DynamicDNAKmer{UInt32})
14

julia> capacity(DynamicAAKmer{UInt8})
0

julia> capacity(DynamicRNAKmer) # NB: UnionAll type
ERROR: MethodError: no method matching capacity(::Type{DynamicRNAKmer})
[...]
```
"""
function capacity(T::Type{<:DynamicKmer{A, U}}) where {A, U}
    bps = BioSequences.bits_per_symbol(A())
    if iszero(bps)
        # If no bits are coding, all bits are used for length.
        # We truncate at typemax(Int)
        clamp(typemax(U), Int)
    else
        div(max_coding_bits(T), bps)
    end
end

BioSequences.encoded_data_eltype(::Type{DynamicKmer{A, U}}) where {A, U} = U

function BioSequences.extract_encoded_element(x::DynamicKmer, i::Integer)
    bps = BioSequences.bits_per_symbol(x)
    shift = 8 * sizeof(x) - (i * bps)
    u = right_shift(x.x, shift)
    mask = one(x.x) << (bps) - one(x.x)
    return u & mask
end

Base.length(x::DynamicKmer) = (x.x & length_mask(typeof(x))) % Int
Base.isempty(x::DynamicKmer) = iszero(x.x)

function Kmer{A, K}(x::DynamicKmer{A}) where {A <: Alphabet, K}
    return @inline derive_type(Kmer{A, K})(x)
end

# This kmer construction and the one below is efficient, as kmers and dynamic kmers
# share a very similar encoding scheme.
function Kmer{A, K, N}(x::DynamicKmer{A}) where {A <: Alphabet, K, N}
    check_kmer(Kmer{A, K, N})
    length(x) == K || error("Must construct kmer from length K DynamicKmer")
    return from_integer(Kmer{A, K, N}, as_integer(x))
end

function Kmer{A1, K, N}(
        x::DynamicKmer{A2}
    ) where {M, A1 <: NucleicAcidAlphabet{M}, A2 <: NucleicAcidAlphabet{M}, K, N}
    check_kmer(Kmer{A1, K, N})
    length(x) == K || error("Must construct kmer from length K DynamicKmer")
    return from_integer(Kmer{A1, K, N}, as_integer(x))
end

const HASH_MASK = 0x6ff6e9f0462d5162 % UInt

Base.copy(x::DynamicKmer) = x
Base.hash(x::DynamicKmer, h::UInt64) = hash(x.x, h ⊻ HASH_MASK) % UInt
fx_hash(x::DynamicKmer, h::UInt64) = ((bitrotate(h, 5) ⊻ x.x) % UInt) * FX_CONSTANT

Base.:(==)(a::DynamicKmer{A, U}, b::DynamicKmer{A, U}) where {A, U} = a === b

# This is similar to Base.promote, except we use it internally only in this package
# for dynamic kmers. We use it to compare dynamic kmers with compatible alphabets,
# which may differ in alphabet or encoded data eltype.
Base.@constprop :aggressive Base.@assume_effects :foldable function promote_dynamic(
        a::DynamicKmer{A, U1},
        b::DynamicKmer{A, U2},
    ) where {A <: Alphabet, U1, U2}
    return if U1 == U2
        (a, b)
    elseif sizeof(U1) < sizeof(U2)
        u = widen_to(U2, a)
        a = _new_dynamic_kmer(A, u)
        (a, b)
    else
        u = widen_to(U1, b)
        b = _new_dynamic_kmer(A, u)
        (a, b)
    end
end

# Same as above, but complicated by the fact that they do not share alphabet.
Base.@constprop :aggressive Base.@assume_effects :foldable function promote_dynamic(
        a::DynamicKmer{<:NucleicAcidAlphabet{N}, U1},
        b::DynamicKmer{<:NucleicAcidAlphabet{N}, U2}
    ) where {N, U1, U2}
    return if U1 == U2
        b = _new_dynamic_kmer(typeof(Alphabet(a)), b.x)
        (a, b)
    elseif sizeof(U1) < sizeof(U2)
        u = widen_to(U2, a)
        a = _new_dynamic_kmer(typeof(Alphabet(b)), u)
        (a, b)
    else
        u = widen_to(U1, b)
        b = _new_dynamic_kmer(typeof(Alphabet(a)), u)
        (a, b)
    end
end

function Base.:(==)(a::DynamicKmer{A, U1}, b::DynamicKmer{A, U2}) where {A, U1, U2}
    (a, b) = promote_dynamic(a, b)
    return a === b
end

function Base.:(==)(a::DynamicKmer{<:NucleicAcidAlphabet{N}, U1}, b::DynamicKmer{<:NucleicAcidAlphabet{N}, U2}) where {N, U1, U2}
    (a, b) = promote_dynamic(a, b)
    return a === b
end

function Base.isless(a::DynamicKmer, b::DynamicKmer)
    (a, b) = promote_dynamic(a, b)
    return isless(a.x, b.x)
end

function Base.cmp(a::DynamicKmer, b::DynamicKmer)
    (a, b) = promote_dynamic(a, b)
    return cmp(a.x, b.x)
end

@inline function Base.getindex(x::DynamicKmer{A}, idx::AbstractUnitRange{<:Integer}) where {A}
    isempty(idx) && return empty(typeof(x))
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

# This is more efficient than the fallback because RC'ing is cheap
function BioSequences.canonical(x::DynamicKmer{<:NucleicAcidAlphabet})
    rc = reverse_complement(x)
    return x < rc ? x : rc
end

function BioSequences._n_gc(x::DynamicKmer{<:NucleicAcidAlphabet})
    u = x.x & ~length_mask(typeof(x))
    return BioSequences.gc_bitcount(u, Alphabet(x))
end

"""
    as_integer(x::DynamicKmer{A, U})::U

Similar to `as_integer` for kmers, but is guaranteed to return a value of `U`,
and the number of coding bits is known at runtime.
"""
function Kmers.as_integer(x::DynamicKmer)
    shift = (8 * sizeof(x) - coding_bits(x))
    return right_shift(x.x, shift)
end

"""
    from_integer(T::Type{<:DynamicKmer{A, U}}, u::U, len::Int)::T

Similar to `from_integer` for `Kmer`, but the length of the resulting `DynamicKmer`
must be passed as an argument. Will error if `len` is larger than the maximal size
supported by `T`.

If `u` is obtained from a `DynamicKmer` with a length different from `len`,
the resulting `DynamicKmer` is reproducible, but not correct and may change between
versions.

# Examples
```jldoctest
julia> d = DynamicDNAKmer{UInt32}(dna"TAGTGCTGTAGGC")
13nt DynamicDNAKmer{UInt32}:
TAGTGCTGTAGGC

julia> u = as_integer(d);

julia> from_integer(typeof(d), u, 13) === d
true

julia> from_integer(typeof(d), u, 12) == d
false
```
"""
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

# Constructor dispatches to RecodingScheme
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
    if iszero(bps)
        _new_dynamic_kmer(typeof(A), len % U)
    end

    if isempty(x)
        return _new_dynamic_kmer(typeof(A), zero(U))
    end
    tup = BioSequences.encoded_data(x)
    # Tuple length has to be at least one, since otherwise kmer would be empty,
    # and branch above would have been taken
    u = if length(tup) == 1
        tup[1] % U
    else
        # Here, we know the tuple has at least 2 elements, so the kmer is at least
        # 128 bits.
        # We checked above that the kmer fits in T, so utype(T) must have at least 128
        # bits, and the utype(T) can contain all bits in the tuple
        # Note also that Kmers.jl only loads on little-endian 64-bit systems,
        # so we can assume those things, too.
        u = zero(utype(T))
        shift = 8 * sizeof(u) - 64
        for i in tup
            u |= left_shift(i % U, shift)
            shift -= 64
        end
    end

    # Kmers keep thier coding bits in the lowest part of the data,
    # and dynamic kmers in the upper.
    # Also, the requested type U may be much bigger than the kmer's tuple,
    # which requires further shift
    shift = 8 * sizeof(U) - len * bps
    u = left_shift(u, shift)

    # Add in length
    u |= len % U
    return _new_dynamic_kmer(typeof(A), u)
end

function build_dynamic_kmer(::Copyable, ::Type{T}, x::DynamicKmer) where {T}
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

# Switch encoding data of `x` to `T`. Error if it doesn't fit.
function switch_backing_encoding(T::Type{<:Unsigned}, x::DynamicKmer{A, U}) where {A, U}
    T == U && return x
    return if sizeof(T) < sizeof(x)
        narrow_to(T, x)
    else
        widen_to(T, x)
    end
end

# Create a DynamicKmer{A, T} containig the same sequence as `x`, efficiently,
# or error if `x` does not fit in that type.
function narrow_to(T::Type{<:Unsigned}, x::DynamicKmer{A, U}) where {A, U}
    newT = DynamicKmer{A, T}
    if max_coding_bits(newT) < length(x)
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

# Create a DynamicKmer{A, T} containig the same sequence as `x`, efficiently
function widen_to(T::Type{<:Unsigned}, x::DynamicKmer{A, U}) where {A, U}
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
4nt DynamicKmer{DNAAlphabet{4}, UInt32}:
AGAY
```
"""
function shift_encoding(x::DynamicKmer{A, U}, encoding::U) where {A <: Alphabet, U <: Unsigned}
    mask = length_mask(typeof(x))
    u = x.x & ~mask
    u = left_shift(u, BioSequences.bits_per_symbol(x))
    u |= left_shift(encoding, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & mask))
end

Base.adjoint(x::DynamicKmer) = x

"""
    @dmer_str -> DynamicKmer

Construct a `DynamicKmer{A, UInt64}` from the given string. The macro must be used with a flag
after the input string, e.g. `d` in `dmer"TAG"d` or `a` in `dmer"PCW"a`, signifying
the alphabet of the dynamic kmer.
The flags `d = DNAAlphabet{2}`, `r = RNAAlphabet{2}` and `a = AminoAcidAlphabet`
are recognized.

# Examples
```jldoctest
julia> dmer"UGCUA"r
5nt DynamicRNAKmer{UInt64}:
UGCUA

julia> dmer"YDLLKKR"a
7aa DynamicAAKmer{UInt64}:
YDLLKKR

julia> dmer"TATTAGCA"d
8nt DynamicDNAKmer{UInt64}:
TATTAGCA
```
"""
macro dmer_str(seq, flag)
    trimmed = BioSequences.remove_newlines(seq)
    # Unlike @dna_str, we default to 2-bit alphabets, because kmers
    # by convention are usually 2-bit only
    return if flag == "dna" || flag == "d"
        DynamicDNAKmer{UInt64}(trimmed)
    elseif flag == "rna" || flag == "r"
        DynamicRNAKmer{UInt64}(trimmed)
    elseif flag == "aa" || flag == "a"
        DynamicAAKmer{UInt64}(trimmed)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

"""
    push(x::T, s)::T where {T <: DynamicKmer}

Create a new `DynamicKmer` of type `T` by adding the symbol `s` to the end of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `ArgumentError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push_first`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DynamicDNAKmer{UInt64}:
TGTGCTGA

julia> d2 = push(d, 'G') # converts from Char to DNA
9nt DynamicDNAKmer{UInt64}:
TGTGCTGAG

julia> d == d2 # does not mutate immutable d
false

julia> push(dmer"RRKRLVD"a, AA_W)
ERROR: ArgumentError: DynamicKmer is already at max capacity
[...]
```
"""
function push(x::DynamicKmer{A, U}, s) where {A, U}
    T = typeof(x)
    E = eltype(x)
    sT = convert(E, s)::E
    enc = U(BioSequences.encode(A(), sT))::U
    bps = BioSequences.bits_per_symbol(A())

    # Update new length. Since length is stored in bottom bits,
    # we can simply add it directly. Neat!
    u = x.x + 0x01
    new_len = (u & length_mask(T)) % Int
    new_len > capacity(T) && throw_argumenterror("DynamicKmer is already at max capacity")

    shift = (8 * sizeof(U)) - (bps * new_len)
    u |= left_shift(enc, shift)
    return _new_dynamic_kmer(A, u)
end

"""
    push_first(x::T, s)::T where {T <: DynamicKmer}

Create a new `DynamicKmer` of type `T` by adding the symbol `s` to the start of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `ArgumentError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DynamicDNAKmer{UInt64}:
TGTGCTGA

julia> d2 = push_first(d, 'G') # converts from Char to DNA
9nt DynamicDNAKmer{UInt64}:
GTGTGCTGA

julia> d == d2 # does not mutate immutable d
false

julia> push_first(dmer"RRKRLVD"a, AA_W)
ERROR: ArgumentError: DynamicKmer is already at max capacity
[...]
```
"""
function push_first(x::DynamicKmer{A, U}, s) where {A, U}
    T = typeof(x)
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
    new_len = (x.x & mask) + 0x01
    (new_len % Int) > capacity(T) && throw_argumenterror("DynamicKmer is already at max capacity")
    u |= new_len
    return _new_dynamic_kmer(A, u)
end

"""
    pop(x::DynamicKmer{A, U})::DynamicKmer{A, U}

Returns a new dynamic kmer with the last symbol of the input `x` removed.
Throws an `ArgumentError` if `x` is empty.

See also: [`pop_first`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"EDEAVY"a
6aa DynamicAAKmer{UInt64}:
EDEAVY

julia> d2 = pop(d)
5aa DynamicAAKmer{UInt64}:
EDEAV

julia> d == d2
false

julia> pop(dmer""a)
ERROR: ArgumentError: Cannot pop empty kmer
[...]
```
"""
function pop(x::DynamicKmer{A, U}) where {A, U}
    isempty(x) && throw_argumenterror("Cannot pop empty kmer")

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
    pop_first(x::DynamicKmer{A, U})::DynamicKmer{A, U}

Returns a new dynamic kmer with the first symbol of the input `x` removed.
Throws an `ArgumentError` if `x` is empty.

See also: [`pop`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"UGCGUAGCUA"r
10nt DynamicRNAKmer{UInt64}:
UGCGUAGCUA

julia> d2 = pop_first(d)
9nt DynamicRNAKmer{UInt64}:
GCGUAGCUA

julia> d == d2
false

julia> pop_first(dmer""r)
ERROR: ArgumentError: Cannot pop empty kmer
[...]
```
"""
function pop_first(x::DynamicKmer{A, U}) where {A, U}
    isempty(x) && throw_argumenterror("Cannot pop empty kmer")

    # Remove length, since we need to shift it to pop first,
    # and shifting would move the length bits
    mask = length_mask(DynamicKmer{A, U})
    u = x.x & ~mask

    # Remove the symbol by shifting
    bps = BioSequences.bits_per_symbol(A())
    u <<= bps

    # Add in length back minus one
    u |= (x.x & mask) - 0x01
    return _new_dynamic_kmer(A, u)
end

@noinline throw_argumenterror(s::String) = throw(ArgumentError(s))
