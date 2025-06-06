"""
    Kmer{A<:Alphabet,K,N} <: BioSequence{A}

An immutable bitstype for representing k-mers - short `BioSequences`
of a fixed length `K`.
Since they can be stored directly in registers, `Kmer`s are generally the most
efficient type of `BioSequence`, when `K` is small and known at compile time.

The `N` parameter is derived from `A` and `K` and is not a free parameter.

See also: [`DNAKmer`](@ref), [`RNAKmer`](@ref), [`AAKmer`](@ref), [`AbstractKmerIterator`](@ref)

# Examples
```jldoctest
julia> RNAKmer{5}("ACGUC")
RNA 5-mer:
ACGUC

julia> Kmer{DNAAlphabet{4}, 6}(dna"TGCTTA")
DNA 6-mer:
TGCTTA

julia> AAKmer{5}((lowercase(i) for i in "KLWYR"))
AminoAcid 5-mer:
KLWYR

julia> RNAKmer{3}("UAUC") # wrong length
ERROR:
[...]
```
"""
struct Kmer{A <: Alphabet, K, N} <: BioSequence{A}
    # The number of UInt is always exactly the number needed, no less, no more.
    # The first symbols pack into the first UInts
    # An UInt with N elements pack into the lowest bits of the UInt, with the
    # first symbols in the higher parts of the UInt.
    # Hence, a sequence A-G of 16-bit elements would pack like:
    # ( ABC, DEFG)
    #  ^ 16 unused bits, the unused bits are always top bits of first UInt
    # Unused bits are always zero

    # This layout complicates some Kmer construction code, but simplifies comparison
    # operators, and we really want Kmers to be efficient.
    data::NTuple{N, UInt}

    # This unsafe method do not clip the head
    function Kmer{A, K, N}(::Unsafe, data::NTuple{N, UInt}) where {A <: Alphabet, K, N}
        check_kmer(Kmer{A, K, N})
        new{A, K, N}(data)
    end
end

# Useful to do e.g. `mer"TAG"d isa Mer{3}`
"""
    Mer{K}

Alias for `Kmer{<:Alphabet, K}`. Useful to dispatch on `K-mers` without regard
for the alphabat

# Example
```jldoctest
julia> mer"DEKR"a isa Mer{4}
true

julia> DNAKmer{6}("TGATCA") isa Mer{6}
true

julia> RNACodon <: Mer{3}
true
```
"""
const Mer{K} = Kmer{<:Alphabet, K}

# Aliases
"Alias for `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K, N} = Kmer{DNAAlphabet{2}, K, N}

"Alias for `Kmer{RNAAlphabet{2},K,N}`"
const RNAKmer{K, N} = Kmer{RNAAlphabet{2}, K, N}

"Alias for `Kmer{AminoAcidAlphabet,K,N}`"
const AAKmer{K, N} = Kmer{AminoAcidAlphabet, K, N}

"Alias for `DNAKmer{3,1}`"
const DNACodon = DNAKmer{3, 1}

"Alias for `RNAKmer{3,1}`"
const RNACodon = RNAKmer{3, 1}

"""
    check_kmer(::Type{Kmer{A,K,N}}) where {A,K,N}

Internal methods that checks that the type parameters are good.

This function should compile to a noop in case the parameterization is good.
"""
@inline function check_kmer(::Type{Kmer{A, K, N}}) where {A, K, N}
    if !(K isa Int)
        throw(ArgumentError("K must be an Int"))
    elseif K < 0
        throw(ArgumentError("Bad kmer parameterisation. K must be greater than 0."))
    end
    n = cld((K * BioSequences.bits_per_symbol(A())) % UInt, (sizeof(UInt) * 8) % UInt) % Int
    if !(N isa Int)
        throw(ArgumentError("N must be an Int"))
    elseif n !== N
        # This has been significantly changed conceptually from before. Now we
        # don't just check K, but *enforce* the most appropriate N for K.
        throw(ArgumentError("Bad kmer parameterisation. For K = $K, N should be $n"))
    end
end

################################################
# Compile-time functions computed on Kmer types
################################################

@inline ksize(::Type{<:Kmer{A, K}}) where {A, K} = K
@inline nsize(::Type{<:Kmer{A, K, N}}) where {A, K, N} = N
@inline n_unused(::Type{<:Kmer{A, K, N}}) where {A, K, N} = capacity(Kmer{A, K, N}) - K
@inline bits_unused(T::Type{<:Kmer}) =
    n_unused(T) * BioSequences.bits_per_symbol(Alphabet(T))

@inline function n_coding_elements(::Type{<:Kmer{A, K}}) where {A, K}
    cld(BioSequences.bits_per_symbol(A()) * K, 8 * sizeof(UInt))
end

@inline function per_word_capacity(::Type{<:Kmer{A}}) where {A}
    div(8 * sizeof(UInt), BioSequences.bits_per_symbol(A()))
end

@inline function capacity(::Type{<:Kmer{A, K, N}}) where {A, K, N}
    per_word_capacity(Kmer{A, K, N}) * N
end

@inline function elements_in_head(::Type{<:Kmer{A, K, N}}) where {A, K, N}
    per_word_capacity(Kmer{A, K, N}) - n_unused(Kmer{A, K, N})
end

"""
    derive_type(::Type{Kmer{A, K}}) -> Type{Kmer{A, K, N}}

Compute the fully parameterized kmer type from only the parameters `A` and `K`.
"""
@inline derive_type(::Type{Kmer{A, K}}) where {A, K} =
    Kmer{A, K, n_coding_elements(Kmer{A, K})}

@inline zero_tuple(T::Type{<:Kmer}) = ntuple(i -> zero(UInt), Val{nsize(T)}())

@inline function zero_kmer(::Type{<:Kmer{A, K}}) where {A, K}
    T2 = derive_type(Kmer{A, K})
    T2(unsafe, zero_tuple(T2))
end

##################
# Various methods
##################

# BioSequences interface
Base.length(x::Kmer) = ksize(typeof(x))
Base.copy(x::Kmer) = x # immutable
BioSequences.encoded_data_eltype(::Type{<:Kmer}) = UInt

# BioSequences helper methods
BioSequences.encoded_data(seq::Kmer) = seq.data

# Misc methods
Base.summary(x::Kmer{A, K, N}) where {A, K, N} = string(eltype(x), ' ', K, "-mer")

function Base.show(io::IO, ::MIME"text/plain", s::Kmer)
    println(io, summary(s), ':')
    print(io, s)
end

# TODO: This is only efficient because the compiler, through Herculean effort,
# is able to completely unroll and inline the indexing operation.
@inline function _cmp(x::Kmer{A1, K1}, y::Kmer{A2, K2}) where {A1, A2, K1, K2}
    if K1 == K2
        cmp(x.data, y.data)
    else
        m = min(K1, K2)
        a = @inline x[1:m]
        b = @inline y[1:m]
        c = cmp(a.data, b.data)
        if iszero(c)
            K1 < K2 ? -1 : K2 < K1 ? 1 : 0
        else
            c
        end
    end
end

# Here, we don't allow comparing twobit to fourbit sequences. We could do this semantically,
# but this would open a whole can of worms, be impossible to optimise and defeat the purpose
# of using Kmers.
Base.cmp(x::Kmer{A}, y::Kmer{A}) where {A} = _cmp(x, y)
Base.cmp(x::Kmer{<:FourBit}, y::Kmer{<:FourBit}) = _cmp(x, y)
Base.cmp(x::Kmer{<:TwoBit}, y::Kmer{<:TwoBit}) = _cmp(x, y)
Base.cmp(x::Kmer{A}, y::Kmer{B}) where {A, B} = throw(MethodError(cmp, (x, y)))

Base.isless(x::Kmer, y::Kmer) = @inline(cmp(x, y)) == -1
Base.:(==)(x::Kmer, y::Kmer) = iszero(@inline cmp(x, y))

Base.:(==)(x::Kmer, y::BioSequence) = throw(MethodError(==, (x, y)))
Base.:(==)(x::BioSequence, y::Kmer) = throw(MethodError(==, (x, y)))

Base.hash(x::Kmer, h::UInt) = hash(x.data, h ⊻ ksize(typeof(x)))

if Sys.WORD_SIZE != 64
    error("Kmer.jl only supports 64-bit systems")
end

# These constants are from the original implementation
@static if Sys.WORD_SIZE == 32
    # typemax(UInt32) / golden ratio
    const FX_CONSTANT = 0x9e3779b9
elseif Sys.WORD_SIZE == 64
    # typemax(UInt64) / pi
    const FX_CONSTANT = 0x517cc1b727220a95
else
    error("Invalid word size")
end

# This implementation is translated from the Rust compiler source code,
# licenced under MIT. The original source is the Firefox source code,
# also freely licensed.
"""
    fx_hash(x, [h::UInt])::UInt

An implementation of `FxHash`. This hash function is extremely fast, but the hashes
are of poor quality compared to Julia's default MurmurHash3. In particular:
* The hash function does not have a good avalanche effect, e.g. the lower bits
  of the result depend only on the top few bits of the input
* The bitpattern zero hashes to zero

However, for many applications, `FxHash` is good enough, if the cost of the
higher rate of hash collisions are offset by the faster speed.

The precise hash value of a given kmer is not guaranteed to be stable across minor
releases of Kmers.jl, but _is_ guaranteed to be stable across minor versions of
Julia.

# Examples
```jldoctest
julia> x = fx_hash(mer"KWQLDE"a);

julia> y = fx_hash(mer"KWQLDE"a, UInt(1));

julia> x isa UInt
true

julia> x == y
false
```
"""
function fx_hash(x::Kmer, h::UInt)
    for i in x.data
        h = (bitrotate(h, 5) ⊻ i) * FX_CONSTANT
    end
    h
end
fx_hash(x) = fx_hash(x, zero(UInt))

Base.adjoint(x::Kmer) = x

"""
    as_integer(x::Kmer)::Unsigned

Get the encoded integer representation of kmer `x`. The returned value is an
unsigned integer between `UInt8` and `UInt128`, with the smallest type chosen
which can fit the coding bits.
Throws an `ArgumentError` if passed kmers with more than 128 bits.

!!! warning
    The value of the encoded representation is an implementation detail, and may
    change in minor versions of Kmers.jl.
    However, the value has the following guarantees:
    * For a given kmer and version of Kmers.jl, the value is deterministic.
    * If the alphabet has unique encodings for each symbol, then two kmers of
      the same length are guaranteed to have distinct encodings.
    * The integer has no more bits than the total number of bits in the encoding
      of the kmer, and the smallest possible unsigned bitstype integer type is
      used.

# Examples
```jldoctest
julia> as_integer(mer"AACT"d)
0x07

julia> as_integer(mer"CT"d) # different K, same encoding
0x07

julia> as_integer(mer"KWPQHVY"a)
0x000b110e05081312

julia> as_integer(mer"VEEKEGVLIKLRK"a)
0x0000001306060b0607130a090b0a010b

julia> Kmers.as_integer(AAKmer{17}("A"^17))
ERROR: ArgumentError: Must have at most 128 bits in encoding
[...]
```
"""
@inline function as_integer(x::Kmer{A, K}) where {A, K}
    isempty(x) && return 0x00
    bits = K * BioSequences.bits_per_symbol(x)
    su = sizeof(UInt)
    t = x.data
    if bits <= 8
        t[1] % UInt8
    elseif bits <= 16
        t[1] % UInt16
    elseif bits <= 32
        t[1] % UInt32
    elseif bits <= 64
        t[1]
    elseif bits <= 128
        # NB: This is NOT equivalent to a reinterpret, because we guarantee that
        # only the lower `bits` bits are set, which does not correspond to the
        # internal representation of a kmer.
        ((t[1] % UInt128) << 64) | (t[2] % UInt128)
    else
        throw(ArgumentError("Must have at most 128 bits in encoding"))
    end
end

"""
    from_integer(T::Type{<:Kmer}, u)::T

Construct a kmer of type `T` from the unsigned bit-integer `u`.
The value of the returned kmer cannot be relied on, but it is guaranteed that
the roundtripping a value produced by `as_integer` will work, and will
result in the same kmer.

`T` must have at most 128 bits coding bits, and `u` must have no more than 128 bits. 
Since not all bits of `u` may be used when constructing the kmer, different
integers may return the same kmer.

# Examples
```jldoctest
julia> kmer = mer"TGATCGTAGAGTGTA"d;

julia> u = as_integer(kmer); typeof(u)
UInt32

julia> from_integer(typeof(kmer), u) === kmer
true
```
"""
function from_integer end

@inline function from_integer(T::Type{<:Kmer{A, K}}, u::Unsigned) where {A, K}
    from_integer(derive_type(T), u)
end

@inline function from_integer(T::Type{<:Kmer{A, K, N}}, u::BitUnsigned) where {A, K, N}
    check_kmer(T)
    bits = K * BioSequences.bits_per_symbol(A())
    iszero(bits) && return zero_kmer(T)
    su = sizeof(u) * 8
    if su > 128
        throw(ArgumentError("Cannot use integers larger than 128 bits"))
    end
    if bits > 128
        throw(ArgumentError("Kmer type must contain at most 128 bits"))
    end
    if bits <= 64
        u = (u % UInt) & get_mask(T)
        return T(unsafe, (u,))
    else
        # NB: This is NOT equivalent to a reinterpret, because we guarantee that
        # only the lower `bits` bits are set, which does not correspond to the
        # internal representation of a kmer.
        # And also because of the masking to ensure no non-coding bits are set.
        a = ((u >> 64) % UInt) & get_mask(T)
        b = u % UInt
        return T(unsafe, (a, b))
    end
end

"""
    push(kmer::Kmer{A, K}, s)::Kmer{A, K+1}

Create a new kmer which is the concatenation of `kmer` and `s`.
Returns a `K+1`-mer.

!!! warn
    Since the output of this function is a `K+1`-mer, use of this function
    in a loop may result in type-instability.

See also: [`push_first`](@ref), [`pop`](@ref), [`shift`](@ref)

# Examples
```jldoctest
julia> push(mer"UGCUGA"r, RNA_G)
RNA 7-mer:
UGCUGAG

julia> push(mer"W"a, 'E')
AminoAcid 2-mer:
WE
```
"""
function push(kmer::Kmer, s)
    bps = BioSequences.bits_per_symbol(kmer)
    A = Alphabet(kmer)
    newT = derive_type(Kmer{typeof(A), length(kmer) + 1})
    # If no free space in data, add new tuple
    new_data = if bits_unused(typeof(kmer)) < bps
        (zero(UInt), kmer.data...)
    else
        kmer.data
    end
    # leftshift_carry the new encoding in.
    encoding = UInt(BioSequences.encode(A, convert(eltype(kmer), s)))
    (_, new_data) = leftshift_carry(new_data, bps, encoding)
    newT(unsafe, new_data)
end

"""
    shift(kmer::Kmer{A, K}, s)::Kmer{A, K}

Push `symbol` onto the end of `kmer`, and pop the first symbol in `kmer`.
Unlike `push`, this preserves the input type, and is less likely to result in
type instability.

See also: [`shift_first`](@ref), [`push`](@ref)

# Examples
```jldoctest
julia> shift(mer"TACC"d, DNA_A)
DNA 4-mer:
ACCA

julia> shift(mer"WKYMLPIIRS"aa, 'F')
AminoAcid 10-mer:
KYMLPIIRSF
```
"""
function shift(kmer::Kmer{A}, s) where {A}
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    shift_encoding(kmer, encoding)
end

"""
    push_first(kmer::Kmer{A, K}, s)::Kmer{A, K+1}

Create a new kmer which is the concatenation of `s` and `kmer`.
Returns a `K+1`-mer. Similar to [`push`](@ref), but places the new symbol `s`
at the front.

!!! warn
    Since the output of this function is a `K+1`-mer, use of this function
    in a loop may result in type-instability.

See also: [`push`](@ref),  [`pop`](@ref), [`shift`](@ref)

# Examples
```jldoctest
julia> push_first(mer"GCU"r, RNA_G)
RNA 4-mer:
GGCU

julia> push_first(mer"W"a, 'E')
AminoAcid 2-mer:
EW
```
"""
function push_first(kmer::Kmer{A}, s) where {A}
    bps = BioSequences.bits_per_symbol(A())
    newT = derive_type(Kmer{A, length(kmer) + 1})
    # If no free space in data, add new tuple
    new_data = if bits_unused(typeof(kmer)) < bps
        (zero(UInt), kmer.data...)
    else
        kmer.data
    end
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    head = first(new_data) | left_shift(encoding, (elements_in_head(newT) - 1) * bps)
    newT(unsafe, (head, tail(new_data)...))
end

"""
    shift_first(kmer::kmer, symbol)::typeof(kmer)

Push `symbol` onto the start of `kmer`, and pop the last symbol in `kmer`.

See also: [`shift`](@ref), [`push`](@ref)

# Examples
```jldoctest
julia> shift_first(mer"TACC"d, DNA_A)
DNA 4-mer:
ATAC

julia> shift_first(mer"WKYMLPIIRS"aa, 'F')
AminoAcid 10-mer:
FWKYMLPIIR
```
"""
function shift_first(kmer::Kmer{A}, s) where {A}
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    shift_first_encoding(kmer, encoding)
end

function shift_first_encoding(kmer::Kmer{A}, encoding::UInt) where {A}
    isempty(kmer) && return kmer
    bps = BioSequences.bits_per_symbol(A())
    (_, new_data) = rightshift_carry(kmer.data, bps, zero(UInt))
    head =
        first(new_data) | left_shift(encoding, (elements_in_head(typeof(kmer)) - 1) * bps)
    typeof(kmer)(unsafe, (head, tail(new_data)...))
end

"""
    pop(kmer::Kmer{A, K})::Kmer{A, K-1}

Returns a new kmer with the last symbol of the input `kmer` removed.
Throws an `ArgumentError` if `kmer` is empty.

!!! warn
    Since the output of this function is a `K-1`-mer, use of this function
    in a loop may result in type-instability.

See also: [`pop_first`](@ref), [`push`](@ref), [`shift`](@ref)

# Examples
```jldoctest
julia> pop(mer"TCTGTA"d)
DNA 5-mer:
TCTGT

julia> pop(mer"QPSY"a)
AminoAcid 3-mer:
QPS

julia> pop(mer""a)
ERROR: ArgumentError:
[...]
```
"""
function pop(kmer::Kmer{A}) where {A}
    isempty(kmer) && throw(ArgumentError("Cannot pop 0-mer"))
    bps = BioSequences.bits_per_symbol(A())
    newT = derive_type(Kmer{A, length(kmer) - 1})
    (_, new_data) = rightshift_carry(kmer.data, bps, zero(UInt))
    new_data = if elements_in_head(typeof(kmer)) == 1
        tail(new_data)
    else
        new_data
    end
    newT(unsafe, new_data)
end

"""
    pop_first(kmer::Kmer{A, K})::Kmer{A, K-1}

Returns a new kmer with the first symbol of the input `kmer` removed.
Throws an `ArgumentError` if `kmer` is empty.

!!! warn
    Since the output of this function is a `K-1`-mer, use of this function
    in a loop may result in type-instability.

See also: [`pop`](@ref), [`push`](@ref), [`shift`](@ref)

# Examples
```jldoctest
julia> pop_first(mer"TCTGTA"d)
DNA 5-mer:
CTGTA

julia> pop_first(mer"QPSY"a)
AminoAcid 3-mer:
PSY

julia> pop_first(mer""a)
ERROR: ArgumentError:
[...]
```
"""
function pop_first(kmer::Kmer{A}) where {A}
    isempty(kmer) && throw(ArgumentError("Cannot pop 0-mer"))
    data = if elements_in_head(typeof(kmer)) == 1
        tail(kmer.data)
    else
        bps = BioSequences.bits_per_symbol(A())
        bits_used = 8 * sizeof(UInt) - (bits_unused(typeof(kmer)) + bps)
        mask = left_shift(UInt(1), bits_used) - UInt(1)
        (first(kmer.data) & mask, tail(kmer.data)...)
    end
    newT = derive_type(Kmer{A, length(kmer) - 1})
    newT(unsafe, data)
end

# Get a mask 0x0001111 ... masking away the unused bits of the head element
# in the UInt tuple
@inline function get_mask(T::Type{<:Kmer})
    UInt(1) << (8 * sizeof(UInt) - bits_unused(T)) - 1
end
