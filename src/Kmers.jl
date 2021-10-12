# Kmers.jl
# ========
#
# Module for kmers and short sequence fragments.
#
# This file is a part of the Kmers.jl, a package in the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

__precompile__()

module Kmers

export
    ###
    ### Mers
    ###

    # Type & aliases
    Kmer,
    DNAKmer,
    DNAKmer27,
    DNAKmer31,
    DNAKmer63,
    RNAKmer,
    RNAKmer27,
    RNAKmer31,
    RNAKmer63,
    AAKmer,
    DNACodon,
    RNACodon,

    # Iteration
    fw_neighbors,
    bw_neighbors,
    
    ###
    ### Sequence literals
    ###
    
    @mer_str,
    @bigmer_str

using BioSequences

include("kmer.jl")

end # module
