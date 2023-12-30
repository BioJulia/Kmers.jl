# Kmers.jl
# ========
#
# Module for kmers and short sequence fragments.
#
# This file is a part of the Kmers.jl, a package in the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/Kmers.jl/blob/master/LICENSE
module Kmers

export Kmer,
    Mer,
    DNAKmer,
    RNAKmer,
    AAKmer,
    DNACodon,
    RNACodon,
    ReverseGeneticCode,
    reverse_translate,
    reverse_translate!,
    @mer_str,

    # Immutable operations
    push,
    push_first,
    shift,
    shift_first,
    pop,
    pop_first,

    # Iterators
    FwKmers,
    FwDNAMers,
    FwRNAMers,
    FwAAMers,
    CanonicalKmers,
    CanonicalDNAMers,
    CanonicalRNAMers,
    UnambiguousKmers,
    UnambiguousDNAMers,
    UnambiguousRNAMers,

    # Reverse translation
    CodonSet,
    delete, # push already exported

    ##################
    # Re-exports
    ##################
    # BioSymbols re-exports
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
    complement,
    reverse_complement,
    canonical,
    iscanonical

# Kmers.jl is tightly coupled to BioSequences and relies on much of its internals.
# Hence, we do not care about carefully importing specific symbols
using BioSequences

# This is a documented method, not internals
using Base: tail

"""
    Kmers.Unsafe

Internal trait object used to access unsafe methods of functions.
`unsafe` is the singleton of `Unsafe`.
"""
struct Unsafe end
const unsafe = Unsafe()

const FourBit = Union{DNAAlphabet{4}, RNAAlphabet{4}}
const TwoBit = Union{DNAAlphabet{2}, RNAAlphabet{2}}
const Bytes = Union{String, SubString{String}, AbstractVector{UInt8}}
const BitInteger =
    Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Int128, UInt128}

include("tuple_bitflipping.jl")
include("kmer.jl")
include("construction.jl")
include("indexing.jl")
include("transformations.jl")
include("revtrans.jl")

include("iterators/common.jl")
include("iterators/FwKmers.jl")
include("iterators/CanonicalKmers.jl")
include("iterators/UnambiguousKmers.jl")
#include("iterators/SpacedKmers.jl")

end # module
