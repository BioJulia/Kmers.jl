###
### Mer specific specializations of src/biosequence/counting.jl
###

for i in [(:_count_a, :a_bitcount), (:_count_c, :c_bitcount), (:_count_g, :g_bitcount), (:_count_t, :t_bitcount)]
    @eval begin
        @inline function $(i[1])(alph::A, head::UInt64, tail...) where {A<:NucleicAcidAlphabet}
            return $(i[2])(head, alph) + $(i[1])(alph, tail...)
        end
        @inline $(i[1])(alph::A) where {A<:NucleicAcidAlphabet} = 0
    end
end

@inline function _count_gc(alph::A, head::UInt64, tail...) where {A<:NucleicAcidAlphabet}
    return gc_bitcount(head, alph) + _count_gc(alph, tail...)
end
@inline _count_gc(::A) where {A<:NucleicAcidAlphabet} = 0

count_a(x::Kmer{A,K,N})  where {A<:NucleicAcidAlphabet,K,N} = _count_a(A(), x.data...) - n_unused(x)
count_c(x::Kmer{A,K,N})  where {A<:NucleicAcidAlphabet,K,N} = _count_c(A(), x.data...)
count_g(x::Kmer{A,K,N})  where {A<:NucleicAcidAlphabet,K,N} = _count_g(A(), x.data...)
count_t(x::Kmer{A,K,N})  where {A<:NucleicAcidAlphabet,K,N} = _count_t(A(), x.data...)

count_gc(x::Kmer{A,K,N}) where {A<:NucleicAcidAlphabet,K,N} = _count_gc(A(), x.data...)
Base.count(::typeof(isGC), x::Kmer{A,K,N}) where {A<:NucleicAcidAlphabet,K,N} = count_gc(x)

# TODO: Expand to Amino Acid Kmers as well...
@inline function Base.count(::typeof(!=), a::Kmer{A,K,N}, b::Kmer{A,K,N}) where {A<:NucleicAcidAlphabet,K,N}
    ad = a.data
    bd = b.data
    sum = 0
    @inbounds for i in 1:N
        sum += BioSequences.mismatch_bitcount(ad[i], bd[i], A())
    end
    return sum
end

# TODO: Expand to Amino Acid Kmers as well...
@inline function Base.count(::typeof(==), a::Kmer{A,K,N}, b::Kmer{A,K,N}) where {A<:NucleicAcidAlphabet,K,N}
    ad = a.data
    bd = b.data
    sum = 0
    @inbounds for i in 1:N
        sum += BioSequences.match_bitcount(ad[i], bd[i], A())
    end
    return sum - n_unused(a)
end