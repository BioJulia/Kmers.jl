"""
    Oligomer{A <: Alphabet, U <: Unsigned} <: BioSequence{A}

Dynamic kmers are immutable, bitstype `BioSequence`s similar to `Kmer`s.
However, unlike the `Kmer` type, the length of a dynamic kmer is a run time
value, and not a compile time value.

Dynamic kmers types have a maximum number of symbols they can store,
see [`capacity`](@ref) for details.

Dynamic kmers are slightly less efficient than regular kmers.
They are useful when a workload includes kmers of varying sizes, where the
length specialization of the `Kmer` type would cause excessive compilation
and/or type instability.

See also: [`DNAOligomer`](@ref), [`Kmer`](@ref)

# Examples
```jldoctest
julia> m = RNAOligomer{UInt32}(rna"AUGCUGA")
7nt RNAOligomer{UInt32}:
AUGCUGA

julia> reverse_complement(m)
7nt RNAOligomer{UInt32}:
UCAGCAU

julia> DNAKmer{7}(m)
DNA 7-mer:
ATGCTGA
```
"""
struct Oligomer{A <: Alphabet, U <: Unsigned} <: BioSequence{A}
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

Base.summary(x::Oligomer{<:Union{DNAAlphabet, RNAAlphabet}}) = string(length(x), "nt ", typeof(x))
Base.summary(x::Oligomer{AminoAcidAlphabet}) = string(length(x), "aa ", typeof(x))
Base.summary(x::Oligomer) = string(length(x), "-symbol ", typeof(x))

function Base.show(io::IO, ::MIME"text/plain", s::Oligomer)
    println(io, summary(s), ':')
    return print(io, s)
end

Base.empty(::Type{<:Oligomer{A, U}}) where {A, U} = _new_dynamic_kmer(A, zero(U))

utype(::Type{<:Oligomer{A, U}}) where {A, U} = U

"Alias for Oligomer{DNAAlphabet{2}, <:Unsigned}"
const DNAOligomer{U} = (Oligomer{DNAAlphabet{2}, U} where {U <: Unsigned})

"Alias for Oligomer{RNAAlphabet{2}, <:Unsigned}"
const RNAOligomer{U} = (Oligomer{RNAAlphabet{2}, U} where {U <: Unsigned})

"Alias for Oligomer{AminoAcidAlphabet, <:Unsigned}"
const AAOligomer{U} = (Oligomer{AminoAcidAlphabet, U} where {U <: Unsigned})

Base.@constprop :aggressive Base.@assume_effects :foldable function max_coding_bits(
        ::Type{Oligomer{A, U}}
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

@inline function length_mask(T::Type{<:Oligomer})
    U = BioSequences.encoded_data_eltype(T)
    return one(U) << (8 * sizeof(U) - max_coding_bits(T)) - one(U)
end

@inline function top_mask(::Type{U}, len::Integer) where {U <: Unsigned}
    return left_shift(typemax(U), (8 * sizeof(U) - len))
end

@inline function top_mask(T::Type{<:Oligomer}, len::Integer)
    return top_mask(BioSequences.encoded_data_eltype(T), len)
end

@inline function coding_bits(x::Oligomer)
    return BioSequences.bits_per_symbol(x) * length(x)
end

@inline function noncoding_bits(x::Oligomer)
    return 8 * sizeof(x) - coding_bits(x)
end

@inline function coding_mask(x::Oligomer)
    return iszero(x.x) ? x.x : top_mask(typeof(x), coding_bits(x))
end

"""
    capacity(T::Type{<:Oligomer{A, U}})::Int

Compute the maximum number of symbols that an instance of the concrete
type `T` can contain.
This computation is a compile time constant and so should not take
any runtime computation.

The value of this number is not guaranteed to be stable across versions of Kmers.

If `B` is the number of bits per symbol of `A`, the answer is `clamp(typemax(U), Int)`,
else the answer is a number in `0:div(8 * sizeof(U), B)`

# Examples
```jldoctest
julia> capacity(DNAOligomer{UInt32})
14

julia> capacity(AAOligomer{UInt8})
0

julia> capacity(RNAOligomer) # NB: UnionAll type
ERROR: MethodError: no method matching capacity(::Type{RNAOligomer})
[...]
```
"""
Base.@constprop :aggressive Base.@assume_effects :foldable function capacity(
        T::Type{<:Oligomer{A, U}}
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

BioSequences.encoded_data_eltype(::Type{Oligomer{A, U}}) where {A, U} = U

function BioSequences.extract_encoded_element(x::Oligomer, i::Integer)
    bps = BioSequences.bits_per_symbol(x)
    shift = 8 * sizeof(x) - (i * bps)
    u = right_shift(x.x, shift)
    mask = one(x.x) << (bps) - one(x.x)
    return u & mask
end

Base.length(x::Oligomer) = (x.x & length_mask(typeof(x))) % Int
Base.isempty(x::Oligomer) = iszero(x.x)

function Kmer{A, K}(x::Oligomer{A}) where {A <: Alphabet, K}
    return @inline derive_type(Kmer{A, K})(x)
end

# This kmer construction and the one below is efficient, as kmers and dynamic kmers
# share a very similar encoding scheme.
function Kmer{A, K, N}(x::Oligomer{A}) where {A <: Alphabet, K, N}
    check_kmer(Kmer{A, K, N})
    length(x) == K || error("Must construct kmer from length K Oligomer")
    return from_integer(Kmer{A, K, N}, as_integer(x))
end

function Kmer{A1, K, N}(
        x::Oligomer{A2}
    ) where {M, A1 <: NucleicAcidAlphabet{M}, A2 <: NucleicAcidAlphabet{M}, K, N}
    check_kmer(Kmer{A1, K, N})
    length(x) == K || error("Must construct kmer from length K Oligomer")
    return from_integer(Kmer{A1, K, N}, as_integer(x))
end

const HASH_MASK = 0x6ff6e9f0462d5162 % UInt

Base.copy(x::Oligomer) = x
Base.hash(x::Oligomer, h::UInt64) = hash(x.x, h ⊻ HASH_MASK) % UInt
fx_hash(x::Oligomer, h::UInt64) = ((bitrotate(h, 5) ⊻ x.x) % UInt) * FX_CONSTANT

Base.:(==)(a::Oligomer{A, U}, b::Oligomer{A, U}) where {A, U} = a === b

# This is similar to Base.promote, except we use it internally only in this package
# for dynamic kmers. We use it to compare dynamic kmers with compatible alphabets,
# which may differ in alphabet or encoded data eltype.
Base.@constprop :aggressive Base.@assume_effects :foldable function promote_dynamic(
        a::Oligomer{A, U1},
        b::Oligomer{A, U2},
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
        a::Oligomer{<:NucleicAcidAlphabet{N}, U1},
        b::Oligomer{<:NucleicAcidAlphabet{N}, U2}
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

function Base.:(==)(a::Oligomer{A, U1}, b::Oligomer{A, U2}) where {A, U1, U2}
    (a, b) = promote_dynamic(a, b)
    return a === b
end

function Base.:(==)(a::Oligomer{<:NucleicAcidAlphabet{N}, U1}, b::Oligomer{<:NucleicAcidAlphabet{N}, U2}) where {N, U1, U2}
    (a, b) = promote_dynamic(a, b)
    return a === b
end

function Base.isless(a::Oligomer, b::Oligomer)
    (a, b) = promote_dynamic(a, b)
    return isless(a.x, b.x)
end

function Base.cmp(a::Oligomer, b::Oligomer)
    (a, b) = promote_dynamic(a, b)
    return cmp(a.x, b.x)
end

@inline function Base.getindex(x::Oligomer{A}, idx::AbstractUnitRange{<:Integer}) where {A}
    isempty(idx) && return empty(typeof(x))
    @boundscheck checkbounds(x, idx)
    bps = BioSequences.bits_per_symbol(x)
    len = length(idx)
    u = left_shift(x.x, (first(idx) - 1) * bps)
    U = BioSequences.encoded_data_eltype(typeof(x))
    u &= top_mask(U, len * bps)
    return _new_dynamic_kmer(A, u | (len % U))
end

function BioSequences.complement(x::Oligomer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    A = typeof(Alphabet(x))
    return _new_dynamic_kmer(A, x.x ⊻ coding_mask(x))
end

function BioSequences.complement(x::Oligomer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}})
    A = typeof(Alphabet(x))
    u = BioSequences.complement_bitpar(x.x, A())
    u &= coding_mask(x)
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

function Base.reverse(x::Oligomer{A}) where {A}
    Bps = BioSequences.BitsPerSymbol(A())
    u = BioSequences.reversebits(x.x, Bps)
    u = left_shift(u, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & length_mask(typeof(x))))
end

function BioSequences.reverse_complement(x::Oligomer{<:NucleicAcidAlphabet})
    return reverse(complement(x))
end

BioSequences.iscanonical(x::Oligomer) = x <= reverse_complement(x)

# This is more efficient than the fallback because RC'ing is cheap
function BioSequences.canonical(x::Oligomer{<:NucleicAcidAlphabet})
    rc = reverse_complement(x)
    return x < rc ? x : rc
end

function BioSequences._n_gc(x::Oligomer{<:NucleicAcidAlphabet})
    u = x.x & ~length_mask(typeof(x))
    return BioSequences.gc_bitcount(u, Alphabet(x))
end

"""
    as_integer(x::Oligomer{A, U})::U

Similar to `as_integer` for kmers, but is guaranteed to return a value of `U`,
and the number of coding bits is known at runtime.
"""
function Kmers.as_integer(x::Oligomer)
    shift = (8 * sizeof(x) - coding_bits(x))
    return right_shift(x.x, shift)
end

"""
    from_integer(T::Type{<:Oligomer{A, U}}, u::U, len::Int)::T

Similar to `from_integer` for `Kmer`, but the length of the resulting `Oligomer`
must be passed as an argument. Will error if `len` is larger than the maximal size
supported by `T`.

If `u` is obtained from a `Oligomer` with a length different from `len`,
the resulting `Oligomer` is reproducible, but not correct and may change between
versions.

# Examples
```jldoctest
julia> d = DNAOligomer{UInt32}(dna"TAGTGCTGTAGGC")
13nt DNAOligomer{UInt32}:
TAGTGCTGTAGGC

julia> u = as_integer(d);

julia> from_integer(typeof(d), u, 13) === d
true

julia> from_integer(typeof(d), u, 12) == d
false
```
"""
function Kmers.from_integer(
        T::Type{Oligomer{A, U}}, x::U, len::Int
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
function Oligomer{T1, U}(x::Oligomer{T2, U}) where {
        B,
        T1 <: NucleicAcidAlphabet{B},
        T2 <: NucleicAcidAlphabet{B},
        U <: Unsigned,
    }
    return _new_dynamic_kmer(T1, x.x)
end

function Oligomer{T1}(x::Oligomer{T2}) where {
        B,
        T1 <: NucleicAcidAlphabet{B},
        T2 <: NucleicAcidAlphabet{B},
    }
    return _new_dynamic_kmer(T1, x.x)
end

# Constructor dispatches to RecodingScheme
function Oligomer{A, U}(x) where {A <: Alphabet, U <: Unsigned}
    return build_dynamic_kmer(RecodingScheme(A(), typeof(x)), Oligomer{A, U}, x)
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
        u
    end

    # Kmers keep thier coding bits in the lowest part of the data,
    # and dynamic kmers in the upper.
    u = left_shift(u, bits_unused(typeof(x)))

    # Add in length
    u |= len % U
    return _new_dynamic_kmer(typeof(A), u)
end

function build_dynamic_kmer(::Copyable, ::Type{T}, x::Oligomer) where {T}
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
function switch_backing_encoding(T::Type{<:Unsigned}, x::Oligomer{A, U}) where {A, U}
    T == U && return x
    return if sizeof(T) < sizeof(x)
        narrow_to(T, x)
    else
        widen_to(T, x)
    end
end

# Create a Oligomer{A, T} containig the same sequence as `x`, efficiently,
# or error if `x` does not fit in that type.
function narrow_to(T::Type{<:Unsigned}, x::Oligomer{A, U}) where {A, U}
    newT = Oligomer{A, T}
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

# Create a Oligomer{A, T} containig the same sequence as `x`, efficiently
function widen_to(T::Type{<:Unsigned}, x::Oligomer{A, U}) where {A, U}
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
    shift_encoding(x::Oligomer{A, U}, encoding::U)::typeof(x)

Add `encoding`, a valid encoding in the alphabet of the `x`,
and of the same integer type as that used in `x`,
to the end of dynamic kmer `x` and discarding the first symbol in `x`.

It is the user's responsibility to ensure that `encoding` is valid.

# Examples
```jldoctest
julia> enc = UInt32(0x0a); # encoding of DNA_Y in 4-bit alphabets

julia> kmer = Oligomer{DNAAlphabet{4}, UInt32}("TAGA");

julia> Kmers.shift_encoding(kmer, enc)
4nt Oligomer{DNAAlphabet{4}, UInt32}:
AGAY
```
"""
function shift_encoding(x::Oligomer{A, U}, encoding::U) where {A <: Alphabet, U <: Unsigned}
    mask = length_mask(typeof(x))
    u = x.x & ~mask
    u = left_shift(u, BioSequences.bits_per_symbol(x))
    u |= left_shift(encoding, noncoding_bits(x))
    return _new_dynamic_kmer(A, u | (x.x & mask))
end

Base.adjoint(x::Oligomer) = x

"""
    @dmer_str -> Oligomer

Construct a `Oligomer{A, UInt64}` from the given string. The macro must be used with a flag
after the input string, e.g. `d` in `dmer"TAG"d` or `a` in `dmer"PCW"a`, signifying
the alphabet of the dynamic kmer.
The flags `d = DNAAlphabet{2}`, `r = RNAAlphabet{2}` and `a = AminoAcidAlphabet`
are recognized.

# Examples
```jldoctest
julia> dmer"UGCUA"r
5nt RNAOligomer{UInt64}:
UGCUA

julia> dmer"YDLLKKR"a
7aa AAOligomer{UInt64}:
YDLLKKR

julia> dmer"TATTAGCA"d
8nt DNAOligomer{UInt64}:
TATTAGCA
```
"""
macro dmer_str(seq, flag)
    trimmed = BioSequences.remove_newlines(seq)
    # Unlike @dna_str, we default to 2-bit alphabets, because kmers
    # by convention are usually 2-bit only
    return if flag == "dna" || flag == "d"
        DNAOligomer{UInt64}(trimmed)
    elseif flag == "rna" || flag == "r"
        RNAOligomer{UInt64}(trimmed)
    elseif flag == "aa" || flag == "a"
        AAOligomer{UInt64}(trimmed)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

"""
    push(x::T, s)::T where {T <: Oligomer}

Create a new `Oligomer` of type `T` by adding the symbol `s` to the end of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `BoundsError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push_first`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DNAOligomer{UInt64}:
TGTGCTGA

julia> d2 = push(d, 'G') # converts from Char to DNA
9nt DNAOligomer{UInt64}:
TGTGCTGAG

julia> d == d2 # does not mutate immutable d
false

julia> push(dmer"RRKRLVD"a, AA_W)
ERROR: BoundsError: attempt to access AAOligomer{UInt64} at index [8]
[...]
```
"""
function push(x::Oligomer{A, U}, s) where {A, U}
    T = typeof(x)

    # Update new length. Since length is stored in bottom bits,
    # we can simply add it directly. Neat!
    u = x.x + 0x01
    new_len = (u & length_mask(T)) % Int

    @boundscheck if new_len > capacity(T)
        boundserror(x, capacity(T) + 1)
    end

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
    push_first(x::T, s)::T where {T <: Oligomer}

Create a new `Oligomer` of type `T` by adding the symbol `s` to the start of `x`.
The argument `s` is converted to the element type of `x` first, so e.g. pushing DNA
to an RNA kmer may work.

Throw an `BoundsError` if `x` is already at max capacity.
See [`capacity`](@ref) to obtain the maximum capacity of `T`.

See also: [`push`](@ref), [`pop`](@ref), [`pop_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"TGTGCTGA"d
8nt DNAOligomer{UInt64}:
TGTGCTGA

julia> d2 = push_first(d, 'G') # converts from Char to DNA
9nt DNAOligomer{UInt64}:
GTGTGCTGA

julia> d == d2 # does not mutate immutable d
false

julia> push_first(dmer"RRKRLVD"a, AA_W)
ERROR: BoundsError: attempt to access AAOligomer{UInt64} at index [8]
[...]
```
"""
function push_first(x::Oligomer{A, U}, s) where {A, U}
    T = typeof(x)

    # Update new length. Since length is stored in bottom bits,
    # we can simply add it directly. Neat!
    u = x.x + 0x01
    new_len = (u & length_mask(T)) % Int

    @boundscheck if new_len > capacity(T)
        boundserror(x, capacity(T) + 1)
    end


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
    pop(x::Oligomer{A, U})::Oligomer{A, U}

Returns a new dynamic kmer with the last symbol of the input `x` removed.
Throws an `BoundsError` if `x` is empty.

See also: [`pop_first`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"EDEAVY"a
6aa AAOligomer{UInt64}:
EDEAVY

julia> d2 = pop(d)
5aa AAOligomer{UInt64}:
EDEAV

julia> d == d2
false

julia> pop(dmer""a)
ERROR: BoundsError: attempt to access AAOligomer{UInt64} at index [0]
[...]
```
"""
function pop(x::Oligomer{A, U}) where {A, U}
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
    pop_first(x::Oligomer{A, U})::Oligomer{A, U}

Returns a new dynamic kmer with the first symbol of the input `x` removed.
Throws an `BoundsError` if `x` is empty.

See also: [`pop`](@ref), [`push`](@ref), [`push_first`](@ref)

# Examples
```jldoctest
julia> d = dmer"UGCGUAGCUA"r
10nt RNAOligomer{UInt64}:
UGCGUAGCUA

julia> d2 = pop_first(d)
9nt RNAOligomer{UInt64}:
GCGUAGCUA

julia> d == d2
false

julia> pop_first(dmer""r)
ERROR: BoundsError: attempt to access RNAOligomer{UInt64} at index [0]
[...]
```
"""
function pop_first(x::Oligomer{A, U}) where {A, U}
    isempty(x) && boundserror(x, 0)

    # Remove length, since we need to shift it to pop first,
    # and shifting would move the length bits
    mask = length_mask(Oligomer{A, U})
    u = x.x & ~mask

    # Remove the symbol by shifting
    bps = BioSequences.bits_per_symbol(A())
    u <<= bps

    # Add in length back minus one
    u |= (x.x & mask) - 0x01
    return _new_dynamic_kmer(A, u)
end

@noinline throw_argumenterror(s::String) = throw(ArgumentError(s))

function Base.setindex(kmer::Oligomer{A, U}, v, i::Integer) where {A, U}
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
        seq::Oligomer{<:Union{DNAAlphabet, RNAAlphabet}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true,
        alternative_start::Bool = false,
    )::AAOligomer

Translate a nucleotide `Oligomer` to a `AAOligomer`.
The type of the result is the smallest `AAOligomer`, which is statically known to
have a capacity large enough to hold the result.

If the result does not fit in the largest known `AAOligomer`, throw an exception.
You can increase the largest known `AAOligomer` by loading the package
BitIntegers.jl.

The arguments other than `seq` are identical to the method with `LongSequence`.

# Examples
```jldoctest
julia> d = Oligomer{DNAAlphabet{4}, UInt64}("TGGCCCGATTGA");

julia> translate(dmer"TGGCCCGATTGA"d)
4aa AAOligomer{UInt128}:
WPD*

julia> translate(DNAOligomer{UInt32}("TGGCCCGATTGA"); alternative_start=true)
4aa AAOligomer{UInt64}:
MPD*
```
"""
function BioSequences.translate(
        seq::Oligomer{<:Union{DNAAlphabet, RNAAlphabet}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true,
        alternative_start::Bool = false,
    )::AAOligomer

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
    if alternative_start
        u |= left_shift(0x0c % U, 8 * sizeof(U) - 8)
    end

    return _new_dynamic_kmer(AminoAcidAlphabet, u)
end

@inline function inbounds_aa_from(
        seq::Oligomer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}},
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
        seq::Oligomer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}},
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
        T::Type{<:Oligomer{<:NucleicAcidAlphabet}}
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
