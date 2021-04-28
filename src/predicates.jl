###
### Mer specific specializations of src/biosequence/predicates.jl
###

Base.cmp(x::T, y::T) where {T<:Kmer} = cmp(x.data, y.data)
Base.:(==)(x::T, y::T) where {T<:Kmer} = x.data == y.data
Base.isless(x::T, y::T) where {T<:Kmer} = isless(x.data, y.data)

# TODO: Ensure this is the right way to go.
# See https://github.com/BioJulia/BioSequences.jl/pull/121#discussion_r475234270
Base.hash(x::Kmer{A,K,N}, h::UInt) where {A,K,N} = hash(x.data, h ⊻ K)