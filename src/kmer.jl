# kmer.jl
# =======
#
# Basic kmer types.
#
# This file is a part of the Kmers.jl, a package in the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract AbstractKmer

for inttype in subtypes(Unsigned)
    nbits = sizeof(inttype) * 8
    typename = "Kmer$nbits"
    @eval begin
        bitstype $nbits $(parse(typename)){T<:NucleicAcid, K} <: AbstractKmer
        typealias $(parse("DNA$typename")){K} $(parse(typename)){DNA, K}
        typealias $(parse("RNA$typename")){K} $(parse(typename)){RNA, K}
    end
end

typealias DNAKmer{K} DNAKmer64{K}
typealias RNAKmer{K} RNAKmer64{K}
typealias DNACodon DNAKmer{3}
typealias RNACodon RNAKmer{3}

# We could also define Kmers as a parametric type such as:
#=
immutable Kmer{T<:NucleicAcid, K, I<:Unsigned} <: Sequence
    data::I
end
=#
