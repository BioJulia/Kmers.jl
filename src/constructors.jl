###
### Constructors for Kmer types
###

#=
These are (hopefully!) very optimised kernel functions for building kmer internal
data from individual elements or from sequences. Kmers themselves are static,
tuple-based structs, and so I really didn't want these functions to create memory
allocations or GC activity through use of vectors an such, for what should be
the creation of a single, rather simple value.
=#

"""
    _build_kmer_data(::Type{Kmer{A,K,N}}, seq::LongSequence{A}, from::Int = 1) where {A,K,N}

Construct a ntuple of the bits data for an instance of a Kmer{A,K,N}.

This particular method is specialised for LongSequences, and for when the Kmer
and LongSequence types used, share the same alphabet, since a lot of encoding /
decoding can be skipped, and the problem is mostly one of shunting bits around.
"""
@inline function _build_kmer_data(::Type{Kmer{A,K,N}}, seq::LongSequence{A}, from::Int = 1) where {A,K,N}
    checkmer(Kmer{A,K,N})
    
    bits_per_sym = BioSequences.bits_per_symbol(A()) # Based on alphabet type, should constant fold.
    n_head = elements_in_head(Kmer{A,K,N}) # Based on kmer type, should constant fold.
    n_per_chunk = per_word_capacity(Kmer{A,K,N}) # Based on kmer type, should constant fold.
    
    if from + K - 1 > length(seq)
        return nothing
    end
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in from:(from + n_head - 1)
        bits = UInt64(BioSequences.extract_encoded_element(seq, i))
        head = (head << bits_per_sym) | bits
    end
    
    # And the rest of the sequence
    idx = Ref(from + n_head)
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            bits = UInt64(BioSequences.extract_encoded_element(seq, idx[]))
            body = (body << bits_per_sym) | bits
            idx[] += 1
        end
        return body
    end
    
    # Put head and tail together
    return (head, tail...)
end

function _construct_from_iterable(::Type{Kmer{A,K,N}}, iter, isize::Base.SizeUnknown) where {A,K,N}
    ## All based on alphabet type of Kmer, so should constant fold.
    checkmer(Kmer{A,K,N})
    ET = eltype(Kmer{A,K,N})
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})

    stateful = Iterators.Stateful(iter)
    head_filled = 0

    # Construct the head.
    head = zero(UInt64)
    next = Ref{Union{Tuple{ET,Nothing},Nothing}}(iterate(stateful))
    while next[] !== nothing && head_filled < n_head
        el, _ = next[]
        sym = convert(ET, el)
        # Encode will throw if it cant encode an element.
        head = (head << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
        head_filled += 1
        next[] = iterate(stateful)
    end

    if head_filled != n_head
        error("Iterator did not yield enough elements to construct $K-mer")
    end

    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        chunk_filled = 0
        body = zero(UInt64)
        while next[] !== nothing && chunk_filled < n_per_chunk
            el, _ = next[]
            sym = convert(ET, el)
            # Encode will throw  if it cant encode an element.
            body = (body << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
            chunk_filled += 1
            next[] = iterate(stateful)
        end

        if chunk_filled != n_per_chunk
            error("Iterator did not yield enough elements to construct $K-mer")
        end

        return body
    end
    
    data = (head, tail...)
    
    return Kmer{A,K,N}(data)
end

function _construct_from_iterable(::Type{Kmer{A,K,N}}, iter, isize::Base.HasLength) where {A,K,N}
    ## All based on alphabet type of Kmer, so should constant fold.
    checkmer(Kmer{A,K,N})
    ET = eltype(Kmer{A,K,N})
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})

    ilen = length(iter)
    if ilen != K
        throw(ArgumentError("seq is not the correct length ($ilen ≠ $K)"))
    end
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in firstindex(iter):(firstindex(iter) + n_head - 1)
        sym = convert(ET, iter[i])
        head = (head << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
    end
    
    # And the rest of the sequence
    idx = Ref(firstindex(iter) + n_head)
    tail = ntuple(Val{N - 1}()) do _
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            sym = convert(ET, iter[i])
            body = (body << bits_per_sym) | UInt64(BioSequences.encode(A(), sym))
            idx[] += 1
        end
        return body
    end
    
    data = (head, tail...)
    
    return Kmer{A,K,N}(data)
end


"""
    Kmer{A,K,N}(itr) where {A,K,N}

Construct a `Kmer{A,K,N}` from an iterable.

The most generic constructor of a Kmer from any iterable type
that implements the bare minimum Iterator interface.

# Examples

```jldoctest
julia> ntseq = LongSequence("TTAGC") # 4-bit DNA alphabet
5nt DNA Sequence:
TTAGC

julia> DNAKmer{5}(ntseq) # 2-Bit DNA alphabet
DNA 5-mer:
TTAGC
```
"""
function Kmer{A,K,N}(itr) where {A,K,N}
    return _construct_from_iterable(Kmer{A,K,N}, itr, Base.IteratorSize(typeof(itr)))
end






"""
    Kmer{A,K,N}(seq::BioSequence{A})

Construct a `Kmer{A,K,N}` from a `BioSequence{A}`.

This particular method is specialised for BioSequences, and for when the Kmer
and BioSequence types used, share the same alphabet, since a lot of encoding /
decoding can be skipped, and the problem is mostly one of shunting bits around.
In the case where the alphabet of the Kmer and the alphabet of the BioSequence
differ, dispatch to the more generic constructor occurs instead.

# Examples

```jldoctest
julia> ntseq = LongSequence{DNAAlphabet{2}}("TTAGC") # 2-bit DNA alphabet
5nt DNA Sequence:
TTAGC

julia> DNAKmer{5}(ntseq) # 2-Bit DNA alphabet
DNA 5-mer:
TTAGC
```
"""
@inline function Kmer{A,K,N}(seq::BioSequence{A}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    
    seqlen = length(seq)
    if seqlen != K
        throw(ArgumentError("seq is not the correct length ($seqlen ≠ $K)"))
    end
    
    ## All based on alphabet type of Kmer, so should constant fold.
    bits_per_sym = BioSequences.bits_per_symbol(A())
    n_head = elements_in_head(Kmer{A,K,N})
    n_per_chunk = per_word_capacity(Kmer{A,K,N})
    
    # Construct the head.
    head = zero(UInt64)
    @inbounds for i in 1:n_head
        bits = UInt64(BioSequences.extract_encoded_element(seq, i))
        head = (head << bits_per_sym) | bits
    end
    
    # And the rest of the sequence
    idx = Ref(n_head + 1)
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for _ in 1:n_per_chunk
            bits = UInt64(BioSequences.extract_encoded_element(seq, idx[]))
            body = (body << bits_per_sym) | bits
            idx[] += 1
        end
        return body
    end
    
    data = (head, tail...)
    
    return Kmer{A,K,N}(data)
end


"""
    Kmer{A,K,n}(x::Kmer{B,K,N}) where {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N,n}

A more optimal method of constructing four-bit encoded nucleic
acid kmers, from two-bit encoded nucleic acid kmers. 
"""
@inline function Kmer{A,K,n}(x::Kmer{B,K,N}) where 
    {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N,n}
    checkmer(Kmer{A,K,n})
    newbits = transcode_bits(x.data, B(), A())
    clippedbits = _drophead(newbits, Val{length(newbits) - n}())
    return Kmer{A,K,n}(clippedbits)
end


"""
    Kmer{A,K}(x::Kmer{B,K,N}) where {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N}

A more optimal method of constructing four-bit encoded nucleic
acid kmers, from two-bit encoded nucleic acid kmers.

Works out appropriate N for the output type as a convenience.
"""
@inline function Kmer{A,K}(x::Kmer{B,K,N}) where
    {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N}
   return kmertype(Kmer{A,K})(x) 
end


"""
    Kmer{A}(x::Kmer{B,K,N}) where {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N}

A more optimal method of constructing four-bit encoded nucleic
acid kmers, from two-bit encoded nucleic acid kmers.

Works out appropriate K and N for output type as a convenience.
"""
@inline function Kmer{A}(x::Kmer{B,K,N}) where
    {A<:NucleicAcidAlphabet{4},B<:NucleicAcidAlphabet{2},K,N}
    return Kmer{A,K}(x)
end








# Convenience version of function above so you don't have to work out correct N.
"""
    Kmer{A,K}(itr) where {A,K}

Construct a `Kmer{A,K,N}` from an iterable.

This is a convenience method which will work out the correct `N` parameter, for
your given choice of `A` & `K`.
"""
@inline function Kmer{A,K}(itr) where {A,K}
    T = kmertype(Kmer{A,K})
    return T(itr)
end

"""
    Kmer{A}(itr) where {A}

Construct a `Kmer{A,K,N}` from an iterable.

This is a convenience method which will work out K from the length of `itr`, and
the correct `N` parameter, for your given choice of `A` & `K`.

!!! warning
    Since this gets K from runtime values, this is gonna be slow!
"""
@inline Kmer{A}(itr) where {A} = Kmer{A,length(itr)}(itr)
@inline Kmer(seq::BioSequence{A}) where A = Kmer{A}(seq)

function Kmer{A1}(seq::BioSequence{A2}) where {A1 <: NucleicAcidAlphabet, A2 <: NucleicAcidAlphabet}
    kmertype(Kmer{A1, length(seq)})(seq)
end

@inline function Kmer{A}(nts::Vararg{Union{DNA, RNA}, K}) where {A <: NucleicAcidAlphabet, K}
    return kmertype(Kmer{A, K})(nts)
end

"""
    Kmer(nts::Vararg{DNA,K}) where {K}

Construct a Kmer from a variable number `K` of DNA nucleotides.

# Examples

```jldoctest
julia> Kmer(DNA_T, DNA_T, DNA_A, DNA_G, DNA_C)
DNA 5-mer:
TTAGC
```
"""
@inline Kmer(nt::DNA, nts::Vararg{DNA}) = DNAKmer((nt, nts...))

"""
    Kmer(nts::Vararg{RNA,K}) where {K}

Construct a Kmer from a variable number `K` of RNA nucleotides.

# Examples

```jldoctest
julia> Kmer(RNA_U, RNA_U, RNA_A, RNA_G, RNA_C)
DNA 5-mer:
UUAGC
```
"""
@inline Kmer(nt::RNA, nts::Vararg{RNA}) = RNAKmer((nt, nts...))


"""
    Kmer(seq::String)

Construct a DNA or RNA kmer from a string.

!!! warning
    As a convenience method, this derives the `K`, `Alphabet`, and `N` parameters
    for the `Kmer{A,K,N}` type from the input string.

# Examples

```jldoctest
julia> Kmer("TTAGC")
DNA 5-mer:
TTAGC
```
"""
@inline function Kmer(seq::String)
    seq′ = BioSequences.remove_newlines(seq)
    hast = false
    hasu = false
    for c in seq′
        hast |= ((c == 'T') | (c == 't'))
        hasu |= ((c == 'U') | (c == 'u'))
    end
    if (hast & hasu) | (!hast & !hasu)
        throw(ArgumentError("Can't detect alphabet type from string"))
    end
    A = ifelse(hast & !hasu, DNAAlphabet{2}, RNAAlphabet{2})
    return Kmer{A,length(seq′)}(seq′)
end


