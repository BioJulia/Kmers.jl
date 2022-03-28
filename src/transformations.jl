
# Bit-parallel element nucleotide complementation
@inline function _complement_bitpar(a::A, head::UInt64, tail...) where {A<:NucleicAcidAlphabet{2}}
    return (~head, _complement_bitpar(a, tail...)...)
end

@inline _complement_bitpar(a::A) where {A<:NucleicAcidAlphabet{2}} = ()

@inline function pushfirst(x::Kmer{A,K,N}, nt) where {A,K,N}
    ntbits = UInt64(BioSequences.encode(A(), nt)) << (62 - (64N - 2K))
    #ntbits = UInt64(@inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 0x01]) << (62 - (64N - 2K))
    return Kmer{A,K,N}(_rightshift_carry(2, ntbits, x.data...))
end

@inline function pushlast(x::Kmer{A,K,N}, nt) where {A,K,N}
    ntbits = UInt64(BioSequences.encode(A(), nt))
    #ntbits = UInt64(@inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 0x01])
    _, newbits = _leftshift_carry(2, ntbits, x.data...)
    return Kmer{A,K,N}(newbits)
end


###
### Transformation methods
###

"""
    complement(seq::T) where {T<:Kmer}

Return a kmer's complement kmer.

# Examples

```jldoctest
julia> complement(Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C))
DNA 5-mer:
AATCG
```
"""
@inline function BioSequences.complement(seq::T) where {T<:Kmer}
    return T(_complement_bitpar(Alphabet(seq), seq.data...))
end

"""
    reverse(seq::Kmer{A,K,N}) where {A,K,N}

Return a kmer that is the reverse of the input kmer.

# Examples

```jldoctest
julia> reverse(Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C))
DNA 5-mer:
CGATT
```
"""
@inline function Base.reverse(seq::Kmer{A,K,N}) where {A,K,N}
#    rdata = _reverse(identity, BioSequences.BitsPerSymbol(seq), seq.data...)
    rdata = _reverse(BioSequences.BitsPerSymbol(seq), seq.data...)
    return Kmer{A,K,N}(rightshift_carry(rdata, 64N - 2K))
end

"""
    reverse_complement(seq::Kmer)

Return the kmer that is the reverse complement of the input kmer.

# Examples

```jldoctest
julia> reverse_complement(Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C))
DNA 5-mer:
GCTAA
```
"""
@inline function BioSequences.reverse_complement(seq::Kmer{A,K,N}) where {A,K,N}
    return complement(reverse(seq))
end

#=
@inline function reverse_complement2(seq::Kmer{A,K,N}) where {A,K,N}
    f = x -> complement_bitpar(x, A())
    rdata = _reverse(f, BioSequences.BitsPerSymbol(seq), seq.data...)
    return Kmer{A,K,N}(rightshift_carry(rdata, 64N - 2K))
end
=#

"""
    canonical(seq::Kmer{A,K,N}) where {A,K,N}

Return the canonical sequence of `seq`.

A canonical sequence is the numerical lesser of a kmer and its reverse complement.
This is useful in hashing/counting sequences in data that is not strand specific,
and thus observing the short sequence is equivalent to observing its reverse complement.

# Examples

```jldoctest
julia> canonical(Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C))
DNA 5-mer:
GCTAA
```
"""
@inline function BioSequences.canonical(seq::Kmer{A,K,N}) where {A,K,N}
    if N < 4
		return min(seq, reverse_complement(seq))
	else
		return iscanonical(seq) ? seq : reverse_complement(seq)
	end
end

###
### Old Mer specific specializations of src/biosequence/transformations.jl
### - not currently transferred to new type.

# TODO: Sort this and decide on transferring to new NTuple based kmers or no.

#=
function swap(x::T, i, j) where {T<:AbstractMer}
    i = 2 * length(x) - 2i
    j = 2 * length(x) - 2j
    b = encoded_data(x)
    x = ((b >> i) ⊻ (b >> j)) & encoded_data_type(x)(0x03)
    return T(b ⊻ ((x << i) | (x << j)))
end



function Random.shuffle(x::T) where {T<:AbstractMer}
    # Fisher-Yates shuffle for mers.
    j = lastindex(x)
    for i in firstindex(x):(j - 1)
        j′ = rand(i:j)
        x = swap(x, i, j′)
    end
    return x
end
=#