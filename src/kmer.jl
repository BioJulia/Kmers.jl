# Notes about Kmers's representation:
# Each element is encoded in the same way as a LongSequence, however the order
# is different. In a Kmer, the elements fill from MSB to LSB, from first to
# last tuple index. Unused bits are always zeroed.
# This layout complicates some Kmer construction code, but simplifies comparison
# operators, and we really want Kmers to be efficient.

"""
    Kmer{A<:Alphabet,K,N} <: BioSequence{A}

A parametric, immutable, bitstype for representing k-mers - short sequences
of a fixed length K.

Since they can be stored directly in registers, `Kmer`s are generally the most
efficient type of `BioSequence`, when `K` is small and known at compile time.
"""
struct Kmer{A <: Alphabet, K, N} <: BioSequence{A}
    # The number of UInt is always exactly the number needed, no less, no more.
    # The first symbols pack into the first UInts
    # An UInt with N elements pack into the lowest bits of the UInt, with the
    # first symbols in the higher parts of the UInt.
    # Hence, a sequence A-G of 16-bit elements would pack like:
    # ( ABC, DEFG)
    #  ^ 16 unused bits, the unused bits are always top bits of first UInt
    data::NTuple{N, UInt}

    # This unsafe method do not clip the head
    function Kmer{A, K, N}(::Unsafe, data::NTuple{N, UInt}) where {A <: Alphabet, K, N}
        check_kmer(Kmer{A, K, N})
        new{A, K, N}(data)
    end

    function Kmer{A, K, N}(data::NTuple{N, UInt}) where {A <: Alphabet, K, N}
        check_kmer(Kmer{A, K, N})
        x = n_unused(Kmer{A, K, N}) * BioSequences.bits_per_symbol(A())
        return new(cliphead(x, data...))
    end
end

# Aliases
"Shortcut for the type `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K, N} = Kmer{DNAAlphabet{2}, K, N}

"Shortcut for the type `Kmer{RNAAlphabet{2},K,N}`"
const RNAKmer{K, N} = Kmer{RNAAlphabet{2}, K, N}

"Shortcut for the type `Kmer{AminoAcidAlphabet,K,N}`"
const AAKmer{K, N} = Kmer{AminoAcidAlphabet, K, N}

"Shorthand for `DNAKmer{3,1}`"
const DNACodon = DNAKmer{3, 1}

"Shorthand for `RNAKmer{3,1}`"
const RNACodon = RNAKmer{3, 1}

"""
    check_kmer(::Type{Kmer{A,K,N}}) where {A,K,N}

Internal method - enforces good kmer type parameterisation.

For a given Kmer{A,K,N} of length K, the number of words used to
represent it (N) should be the minimum needed to contain all K symbols.

This function should compile to a noop in case the parameterization is good.
"""
@inline function check_kmer(::Type{Kmer{A, K, N}}) where {A, K, N}
    if !(K isa Int)
        throw(ArgumentError("K must be an Int"))
    elseif K < 1
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

@inline ksize(::Type{<:Kmer{A, K, N}}) where {A, K, N} = K
@inline nsize(::Type{<:Kmer{A, K, N}}) where {A, K, N} = N
@inline n_unused(::Type{<:Kmer{A, K, N}}) where {A, K, N} = capacity(Kmer{A, K, N}) - K

@inline function n_coding_elements(::Type{<:Kmer{A, K}}) where {A, K}
    cld(BioSequences.bits_per_symbol(A()) * K, 8 * sizeof(UInt))
end

@inline function per_word_capacity(::Type{<:Kmer{A}}) where A
    div(8 * sizeof(UInt), BioSequences.bits_per_symbol(A()))
end

@inline function capacity(::Type{<:Kmer{A, K, N}}) where {A, K, N}
    per_word_capacity(Kmer{A, K, N}) * N
end

@inline function elements_in_head(::Type{<:Kmer{A, K, N}}) where {A, K, N}
    per_word_capacity(Kmer{A, K, N}) - n_unused(Kmer{A, K, N})
end

################################################
# Constructors
################################################

zero_tuple(T::Type{<:Kmer}) = ntuple(i -> zero(UInt), nsize(T))

# Generic, unknown size
@inline function construct_generic(::Base.SizeUnknown, T::Type{<:Kmer{A, K}}, itr) where {A, K}
    check_kmer(T)
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for (i, element) in enumerate(itr)
        i > K && error("Length of sequence must be K elements to build Kmer")
        symbol = convert(eltype(A), element)
        carry = UInt(BioSequences.encode(A(), symbol))
        data = leftshift_carry(data, nbits, carry)
    end
    T(unsafe, data)
end

# Generic, size known
@inline function construct_generic_unchecked(::Union{Base.HasLength, Base.HasShape}, T::Type{<:Kmer{A}}, itr) where A
    check_kmer(T)
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for element in itr
        symbol = convert(eltype(A), element)
        carry = UInt(BioSequences.encode(A(), symbol))
        data = leftshift_carry(data, nbits, carry)
    end
    T(unsafe, data)
end

# Generic, size known but length not checked.
@inline function construct_generic(iT::Union{Base.HasLength, Base.HasShape}, T::Type{<:Kmer{A, K}}, itr) where {A, K}
    length(s) == K || error("Length of sequence must be K elements to build Kmer")
    construct_generic_unchecked(iT, T, itr)
end

# BioSequences with the same Alphabet and these element types do not need to decode
# and encode, but can copy the raw bits directly into the kmer
@inline function construct_unchecked(
    T::Type{<:Kmer{A}}, s::BioSequence{A}, data_eltype::Type{E}
) where {A, E <: Union{UInt8, UInt16, UInt32, UInt}}
    check_kmer(T)
    data = zero_tuple(T)
    nbits = BioSequences.bits_per_symbol(A())
    for i in 1:K
        data = leftshift_carry(data, nbits, BioSequences.extract_encoded_element(s, i) % UInt)
    end
    T(unsafe, data)
end

# BioSequence with another element type fall back to the generic length constructor
@inline function construct_unchecked(T::Type{<:Kmer}, s::BioSequence, data_eltype::Type)
    construct_generic_unchecked(Base.HasLength(), T, s)
end

# BioSequence must implement length so we don't need to dispatch on that.
# However, if the encoded data eltype is an unsigned, we can use a specialized method where we don't
# decode each symbol but simply move the encoded data directly into the tuple
function Kmer{A, K, N}(s::BioSequence) where {A, K, N}
    length(s) == K || error("Length of sequence must be K elements to build Kmer")
    construct_unchecked(T, s, BioSequences.encoded_data_eltype(typeof(s)))
end

# Generic constructor: Dispatch on the iteratorsize
function Kmer{A, K, N}(itr) where {A, K, N}
    construct_generic(Base.IteratorSize(typeof(itr)), Kmer{A, K, N}, itr)
end

# To avoid having the O(N) length check. TODO: Use optimised method
function Kmer{A, K, N}(s::Union{String, SubString{String}}) where {A, K, N}
    construct_generic(Base.SizeUnknown(), Kmer{A, K, N}, s)
end

################################################
# Derived constructors
################################################

# Where the parameters of the kmer is not specified in the constructor
function Kmer(s::BioSequence{A}) where A
    K = length(s)
    N = n_coding_elements(Kmer{A, K})
    Kmer{A, K, N}(s)
end

# Where A, but not K is specified
function Kmer{A}(s::Union{String, SubString{String}}) where A
    K = length(s)
    N = n_coding_elements(Kmer{A, K})
    construct_generic_unchecked(Base.HasLength(), Kmer{A, K, N}, s)
end

# TODO: Constructor from LongSequence and LongSubSeq
# where whole coding elements can be copied directly over
# without extracting individual elements

# TODO: Kmer => LongSequence constructor, same as above but opposite, kinda.

# TODO: Constructor from String that predicts the alphabet?
# Maybe implement the guessparse function in BioSequences.jl
# (See related issue), then call it from here.

################################################
# String literals
################################################

macro mer_str(seq, flag)
    trimmed = BioSequences.remove_newlines(seq)
    # Unlike @dna_str, we default to 2-bit alphabets, because kmers
    # by convention are usually 2-bit only
    if flag == "dna" || flag == "d"
        Kmer{DNAAlphabet{2}}(trimmed)
    elseif flag == "rna" || flag == "r"
        Kmer{RNAAlphabet{2}}(trimmed)
    elseif flag == "aa" || flag == "a"
        Kmer{AminoAcidAlphabet}(trimmed)
    else
        error("Invalid type flag: '$(flag)'")
    end
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

function Base.print(io::IO, s::Kmer)
    # TODO: Can be optimised but whatever
    print(io, LongSequence(s))
end
