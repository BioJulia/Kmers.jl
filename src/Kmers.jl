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

import BioSequences:
    BitsPerSymbol

#import BioSymbols

include("kmer.jl")

end # module
