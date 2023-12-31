"""
    SpacedKmers{A <: Alphabet, K, J, S} <: AbstractKmerIterator{A, K}

Iterator of kmers with step size. `J` signifies the step size, `S`
the type of the underlying sequence, and the eltype of the iterator
is `Kmer{A, K, N}` with the appropriate `N`

See also: [`each_dna_codon`](@ref), [`FwKmers`](@ref)

# Examples:
```jldoctest
julia> collect(SpacedDNAMers{3, 2}("AGCGTATA"))
3-element Vector{Kmer{DNAAlphabet{2}, 3, 1}}:
 AGC
 CGT
 TAT
```
"""
struct SpacedKmers{A <: Alphabet, K, J, S} <: AbstractKmerIterator{A, K}
    seq::S

    function SpacedKmers{A, K, J, S}(seq::S) where {A, K, J, S}
        K isa Int || error("K must be an Int")
        K > 0 || error("K must be at least 1")
        J isa Int || error("J must be an Int")
        J > 0 || error("J must be at least 1")
        new{A, K, J, S}(seq)
    end
end

source_type(::Type{SpacedKmers{A, K, J, S}}) where {A, K, J, S} = S
stepsize(::SpacedKmers{A, K, J}) where {A, K, J} = J

@inline function Base.length(it::SpacedKmers{A, K, J}) where {A, K, J}
    src = used_source(RecodingScheme(A(), source_type(typeof(it))), it.seq)
    L = length(src)
    L < K ? 0 : div((L - K), J) + 1
end

SpacedKmers{A, K, J}(s) where {A <: Alphabet, K, J} = SpacedKmers{A, K, J, typeof(s)}(s)

"`SpacedDNAMers{K, J, S}`: Alias for `SpacedKmers{DNAAlphabet{2}, K, J, S}`"
const SpacedDNAMers{K, J, S} = SpacedKmers{DNAAlphabet{2}, K, J, S}

"`SpacedRNAMers{K, J, S}`: Alias for `SpacedKmers{RNAAlphabet{2}, K, J, S}`"
const SpacedRNAMers{K, J, S} = SpacedKmers{RNAAlphabet{2}, K, J, S}

"`SpacedAAMers{K, J, S}`: Alias for `SpacedKmers{AminoAcidAlphabet, K, J, S}`"
const SpacedAAMers{K, J, S} = SpacedKmers{AminoAcidAlphabet, K, J, S}

# TODO: Do we need two function names for this?...
# Could it be each_codon(DNA, s) <- requires RNA/DNA to be exported
# each(DNACodon) <- awkward since DNACodon is merely an alias
"""
    each_dna_codon(s)

Construct an iterator of `DNACodon` from `s`, iterating over every `DNACodon`
in `s`, in-frame, i.e. with a step size of 3.

See also: [`SpacedKmers`](@ref)

Examples:
```jldoctest
julia> collect(each_dna_codon("TGACGATCGAC"))
3-element Vector{Kmer{DNAAlphabet{2}, 3, 1}}:
 TGA
 CGA
 TCG
```
"""
@inline each_dna_codon(s) = SpacedDNAMers{3, 3}(s)

"The `RNA` equivalent of [`each_dna_codon`](@ref)"
@inline each_rna_codon(s) = SpacedRNAMers{3, 3}(s)

@inline function Base.iterate(it::SpacedKmers{A}, state...) where {A}
    iterate_kmer(RecodingScheme(A(), source_type(typeof(it))), it, state...)
end

# TODO: Maybe in all kmer iterators, instantiate it with the source type,
# so we don't have to get the source type in functions (and thus
# it is allwoed to be a costly operation).
# However, this means we instantiate e.g. a FwKmers{A, K, S} and change S
# in the source type in the constructor
@inline function iterate_kmer(
    R::RecodingScheme,
    it::SpacedKmers{A, K},
) where {A <: Alphabet, K}
    length(it.seq) < K && return nothing
    kmer = unsafe_extract(
        R,
        eltype(it),
        used_source(RecodingScheme(A(), source_type(typeof(it))), it.seq),
        1,
    )
    next_index = 1 + max(stepsize(it), K)
    (kmer, (kmer, next_index))
end

# Here, we need to convert to an abstractvector
# TODO: This function and the one above can be merged with the FwKmers one?
@inline function iterate_kmer(
    R::AsciiEncode,
    it::SpacedKmers{A, K, J, S},
) where {A <: Alphabet, K, J, S <: Bytes}
    src = used_source(RecodingScheme(A(), S), it.seq)
    Base.require_one_based_indexing(src)
    length(src) < K && return nothing
    kmer = unsafe_extract(R, eltype(it), src, 1)
    next_index = 1 + max(stepsize(it), K)
    (kmer, (kmer, next_index))
end

@inline function iterate_kmer(
    ::RecodingScheme,
    it::SpacedKmers{A, K, J, S},
    state,
) where {A, K, S, J}
    src = used_source(RecodingScheme(A(), S), it.seq)
    R = RecodingScheme(A(), S)
    Base.require_one_based_indexing(src)
    (kmer, i) = state
    i > lastindex(src) - min(K, J) + 1 && return nothing
    next_i = i + min(K, J)
    # This branch should be resolved statically
    if J ≥ K
        kmer = unsafe_extract(R, eltype(it), src, i)
    else
        for _ in 1:J
            kmer = update_kmer(R, it, kmer, i)
            i += 1
        end
    end
    (kmer, (kmer, next_i))
end

# TODO: Can this function be used more generically, by the other iterators?
# I.e. this simply fetches a single element
@inline function update_kmer(::GenericRecoding, it::SpacedKmers, kmer::Kmer, i::Int)
    symbol = @inbounds it.seq[i]
    shift(kmer, convert(eltype(kmer), symbol))
end

@inline function update_kmer(::Copyable, it::SpacedKmers, kmer::Kmer, i::Int)
    shift_encoding(kmer, UInt(BioSequences.extract_encoded_element(it.seq, i))::UInt)
end

@inline function update_kmer(::TwoToFour, it::SpacedKmers, kmer::Kmer, i::Int)
    encoding = left_shift(UInt(1), UInt(BioSequences.extract_encoded_element(it.seq, i)))
    shift_encoding(kmer, encoding)
end

@inline function update_kmer(::FourToTwo, it::SpacedKmers, kmer::Kmer, i::Int)
    encoding = UInt(BioSequences.extract_encoded_element(it.seq, i))::UInt
    if count_ones(encoding) != 1
        throw(
            BioSequences.EncodeError(
                Alphabet(kmer),
                reinterpret(eltype(it.seq), encoding % UInt8),
            ),
        )
    end
    shift_encoding(kmer, trailing_zeros(encoding) % UInt)
end

@inline function update_kmer(::AsciiEncode, it::SpacedKmers, kmer::Kmer, i::Int)
    src = used_source(RecodingScheme(Alphabet(eltype(it)), source_type(typeof(it))), it.seq)
    Base.require_one_based_indexing(src)
    byte = @inbounds src[i]
    encoding = BioSequences.ascii_encode(Alphabet(eltype(it)), byte)
    if encoding > 0x7f
        throw(BioSequences.EncodeError(Alphabet(eltype(it)), repr(byte)))
    end
    shift_encoding(kmer, encoding % UInt)
end
