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
    # Types & aliases
    NucleicAcid,
    DNA,
    RNA,
    DNA_A,
    DNA_C,
    DNA_G,
    DNA_T,
    DNA_M,
    DNA_R,
    DNA_W,
    DNA_S,
    DNA_Y,
    DNA_K,
    DNA_V,
    DNA_H,
    DNA_D,
    DNA_B,
    DNA_N,
    DNA_Gap,
    ACGT,
    ACGTN,
    RNA_A,
    RNA_C,
    RNA_G,
    RNA_U,
    RNA_M,
    RNA_R,
    RNA_W,
    RNA_S,
    RNA_Y,
    RNA_K,
    RNA_V,
    RNA_H,
    RNA_D,
    RNA_B,
    RNA_N,
    RNA_Gap,
    ACGU,
    ACGUN,
    AminoAcid,
    AA_A,
    AA_R,
    AA_N,
    AA_D,
    AA_C,
    AA_Q,
    AA_E,
    AA_G,
    AA_H,
    AA_I,
    AA_L,
    AA_K,
    AA_M,
    AA_F,
    AA_P,
    AA_S,
    AA_T,
    AA_W,
    AA_Y,
    AA_V,
    AA_O,
    AA_U,
    AA_B,
    AA_J,
    AA_Z,
    AA_X,
    AA_Term,
    AA_Gap,

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
