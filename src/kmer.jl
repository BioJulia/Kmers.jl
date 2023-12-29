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
RNA 5-mer
ACGUC

julia> Kmer{DNAAlphabet{4}, 6}(dna"TGCTTA")
DNA 6-mer
TGCTTA

julia> AAKmer{5}((lowercase(i) for i in "KLWYR"))
AminoAcid 5-mer
TGCTTA

julia> RNAKmer{3}("UAUC") # wrong length
ERROR:
[ ... ]
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

julia> DNAKmer{2}("TGATCA") isa Mer{6}
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

@inline BioSequences.Alphabet(::Kmer{A}) where {A} = A()

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

@inline derive_type(::Type{Kmer{A, K}}) where {A, K} =
    Kmer{A, K, n_coding_elements(Kmer{A, K})}

zero_tuple(T::Type{<:Kmer}) = ntuple(i -> zero(UInt), Val{nsize(T)}())

# TODO: Should this somehow throw a MethodError if N is already parameterized?
function zero_kmer(T::Type{Kmer{A, K}}) where {A, K}
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

@inline function _cmp(x::Kmer{A1, K1}, y::Kmer{A2, K2}) where {A1, A2, K1, K2}
    if K1 < K2
        -1
    elseif K2 < K1
        1
    else
        cmp(x.data, y.data)
    end
end

# Here, we don't allow comparing twobit to fourbit sequences. We could do this semantically,
# but this would open a whole can of worms, be impossible to optimise and defeat the purpose
# of using Kmers.
Base.cmp(x::Kmer{A}, y::Kmer{A}) where {A} = _cmp(x, y)
Base.cmp(x::Kmer{<:FourBit}, y::Kmer{<:FourBit}) = _cmp(x, y)
Base.cmp(x::Kmer{<:TwoBit}, y::Kmer{<:TwoBit}) = _cmp(x, y)

Base.isless(x::Kmer, y::Kmer) = cmp(x, y) == -1
Base.:(==)(x::Kmer, y::Kmer) = iszero(cmp(x, y))

Base.:(==)(x::Kmer, y::BioSequence) = throw(MethodError(==, (x, y)))
Base.:(==)(x::BioSequence, y::Kmer) = throw(MethodError(==, (x, y)))

Base.hash(x::Kmer{A, K, N}, h::UInt) where {A, K, N} = hash(x.data, h ⊻ K)

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
julia> shift(mer"UGCUGA"r, RNA_G)
RNA 7-mer
UGCUGAG

julia> shift(mer"W"a, 'E')
AminoAcid 2-mer
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
DNA 4-mer
ACCA

julia> shift(mer"WKYMLPIIRS"aa, 'F')
AminoAcid 10-mer
KYMLPIIRSF
```
"""
function shift(kmer::Kmer{A}, s) where {A}
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    shift_encoding(kmer, encoding)
end

@inline function shift_encoding(kmer::Kmer, encoding::UInt)
    bps = BioSequences.bits_per_symbol(kmer)
    (_, new_data) = leftshift_carry(kmer.data, bps, encoding)
    (head, tail...) = new_data
    typeof(kmer)(unsafe, (head & get_mask(typeof(kmer)), tail...))
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
julia> shift(mer"GCU"r, RNA_G)
RNA 4-mer
GGCU

julia> shift(mer"W"a, 'E')
AminoAcid 2-mer
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
    (head, tail...) = new_data
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    head |= left_shift(encoding, (elements_in_head(newT) - 1) * bps)
    newT(unsafe, (head, tail...))
end

"""
    shift_first(kmer::kmer, symbol)::typeof(kmer)

Push `symbol` onto the start of `kmer`, and pop the last symbol in `kmer`.

See also: [`shift`](@ref), [`push`](@ref)

# Examples
```jldoctest
julia> shift_first(mer"TACC"d, DNA_A)
DNA 4-mer
ATAC

julia> shift_first(mer"WKYMLPIIRS"aa, 'F')
AminoAcid 10-mer
FWKYMLPIIR
```
"""
function shift_first(kmer::Kmer{A}, s) where {A}
    encoding = UInt(BioSequences.encode(A(), convert(eltype(kmer), s)))
    shift_first_encoding(kmer, encoding)
end

function shift_first_encoding(kmer::Kmer{A}, encoding::UInt) where {A}
    bps = BioSequences.bits_per_symbol(A())
    (_, new_data) = rightshift_carry(kmer.data, bps, zero(UInt))
    (head, tail...) = new_data
    head |= left_shift(encoding, (elements_in_head(typeof(kmer)) - 1) * bps)
    typeof(kmer)(unsafe, (head, tail...))
end

"""
    pop(kmer::Kmer{A, K})::Kmer{A, K-1}

Returns a new kmer with the last symbol of the input `kmer` removed.
Throws an `ArgumentError` if `kmer` is empty.

!!! warn
    Since the output of this function is a `K+1`-mer, use of this function
    in a loop may result in type-instability.

See also: [`push`](@ref), [`shift`](@ref)

# Examples
```jldoctest
julia> pop(mer"TCTGTA"d)
DNA 5-mer
TCTGT

julia> pop(mer"QPSY"a)
AminoAcid 3-mer
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
        (head, tail...) = new_data
        tail
    else
        new_data
    end
    newT(unsafe, new_data)
end

# Get a mask 0x0001111 ... masking away the unused bits of the head element
# in the UInt tuple
@inline function get_mask(T::Type{<:Kmer})
    UInt(1) << (8 * sizeof(UInt) - bits_unused(T)) - 1
end
