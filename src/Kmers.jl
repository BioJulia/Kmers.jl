# Kmers.jl
# ========
#
# Module for kmers and short sequence fragments.
#
# This file is a part of the Kmers.jl, a package in the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/Kmers.jl/blob/master/LICENSE

__precompile__()

module Kmers

export
    # BioSymbols re-exports.
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
    
    # BioSequences re-exports
    Alphabet,
    BioSequence,
    NucleicAcidAlphabet,
    AminoAcidAlphabet,
    DNAAlphabet,
    RNAAlphabet,
    translate,
    
    
    ###
    ### Mers
    ###

    # Type & aliases
    Kmer,
    DNAKmer,
    DNA27mer,
    DNA31mer,
    DNA63mer,
    RNAKmer,
    RNA27mer,
    RNA31mer,
    RNA63mer,
    AAKmer,
    DNACodon,
    RNACodon,

    # Iteration
    EveryKmer,
    SpacedKmers,
    EveryCanonicalKmer,
    SpacedCanonicalKmers,
    fw_neighbors,
    bw_neighbors,

    # Immutable operators
    push,

    # Translation
    reverse_translate,
    reverse_translate!,
    ReverseGeneticCode,
    
    ###
    ### Sequence literals
    ###
    
    @mer_str,
    @bigmer_str

using BioSequences

ispermitted(::DNAAlphabet{2}, nt::DNA) = count_ones(nt) == 1 && isvalid(nt)
ispermitted(::DNAAlphabet{2}, data::UInt) = data < UInt(4)
ispermitted(::DNAAlphabet{4}, nt::DNA) = isvalid(nt)
ispermitted(::DNAAlphabet{4}, data::UInt) = isvalid(DNA, data)
ispermitted(::AminoAcidAlphabet, aa::AminoAcid) = reinterpret(UInt8, aa) <= reinterpret(UInt8, AA_Gap)
ispermitted(::AminoAcidAlphabet, data::UInt) = data <= 0x1b

include("kmer.jl")

include("kmer_iteration/AbstractKmerIterator.jl")
include("kmer_iteration/EveryKmer.jl")
include("kmer_iteration/SpacedKmers.jl")
include("kmer_iteration/EveryCanonicalKmer.jl")
include("kmer_iteration/SpacedCanonicalKmers.jl")

end # module
