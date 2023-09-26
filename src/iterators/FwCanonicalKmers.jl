struct FwCanonicalKmers{A <: Alphabet, K, S} <: AbstractKmerIterator{A, K}
    it::FwKmers{A, K, S}
end

const SameFwCanonicalKmers{A, K, S} = FwCanonicalKmers{S, A, K} where {A, S <: BioSequence{A}}

function FwCanonicalKmers{K}(s) where K
    S = typeof(s)
    A = typeof(Alphabet(S))
    it = FwKmers{S, A, K}(s)
    FwCanonicalKmers{S, A, K}(it)
end

function FwCanonicalKmers{A, K}(s::S) where {S <: BioSequence, A <: Alphabet, K}
    FwCanonicalKmers{S, A, K}(FwKmers{A, K}(s))
end

function FwCanonicalKmers{A, K}(s::S) where {S <: Union{String, SubString{String}}, A <: Alphabet, K}
    s2 = codeunits(s)
    FwCanonicalKmers{typeof(s2), A, K}(s2)
end

Base.IteratorSize(::Type{<:SameFwCanonicalKmers}) = Base.HasLength()
Base.IteratorSize(::Type{<:FwCanonicalKmers{<:BioSequence{<:TwoBit}, <:FourBit}}) = Base.HasLength()
Base.length(it::SameFwCanonicalKmers) = length(it.it)

# Generic iterator for the first element: I think we can do no better than to reverse-complement
# the entire kmer. However, for the following iterations, it's faster to add a single basepair to
# the RC'd kmer than to RC it from scratch, hence we need specialized methods for efficient RC'ing
# of individual bases.
function Base.iterate(it::FwCanonicalKmers{S, A, K}) where {S, A, K}
    itval = iterate(it.it)
    itval === nothing && return nothing
    fw = first(itval)
    rv = reverse_complement(fw)
    (min(fw, rv), (fw, rv, K+1))
end

# Generic fallback
function Base.iterate(
    it::FwCanonicalKmers{S, A, K},
    state::Tuple{Kmer, Kmer, Integer}
) where {S, A, K}
    seq = it.it.seq
    (fw, rv, i) = state
    i > length(seq) && return nothing
    symbol = @inbounds seq[i]
    encoding = UInt(BioSequences.encode(A, symbol))
    rc_encoding = UInt(BioSequences.encode(A, complement(symbol)))
    fw = shift_encoding(fw, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rv), (fw, rv, i+1))
end

# Special method for 2bit -> 2bit
function Base.iterate(
    it::FwCanonicalKmers{S, A, K},
    state::Tuple{Kmer, Kmer, Integer}
) where {K, A <: TwoBit, S <: BioSequence{A}}
    seq = it.it.seq
    (fw, rv, i) = state
    i > length(seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(seq, i))
    rc_encoding = encoding ⊻ 0x3
    fw = shift_encoding(fw, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rv), (fw, rv, i+1))
end

# Special method for 2bit -> 4bit
function Base.iterate(
    it::FwCanonicalKmers{S, A, K},
    state::Tuple{Kmer, Kmer, Integer}
) where {K, A <: FourBit, S <: BioSequence{A}}
    seq = it.it.seq
    (fw, rv, i) = state
    i > length(seq) && return nothing
    encoding = UInt(BioSequences.extract_encoded_element(seq, i))
    # Reverse-complementing like this is surprisingly inefficient.
    # We may want to consider either using a 16-element LUT, or 
    # else simply changing the algorithm such that the whole kmer
    # is reverse-complemented at every iteration
    rc_encoding = reinterpret(UInt8, complement(reinterpret(DNA, encoding % UInt8))) % UInt
    fw = shift_encoding(fw, encoding)
    rv = shift_first_encoding(rv, rc_encoding)
    (min(fw, rv), (fw, rv, i+1))
end
