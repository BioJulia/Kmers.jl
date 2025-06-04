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
    fx_hash,
    derive_type,
    as_integer,

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
    FwRvIterator,
    CanonicalKmers,
    CanonicalDNAMers,
    CanonicalRNAMers,
    UnambiguousKmers,
    UnambiguousDNAMers,
    UnambiguousRNAMers,
    SpacedKmers,
    SpacedDNAMers,
    SpacedRNAMers,
    SpacedAAMers,
    each_codon,

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
    AminoAcid,

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

# TODO: Remove this ugly hack when 1.11 becomes LTS
if VERSION >= v"1.11.0-DEV.469"
    let str = """
        public unsafe_shift_from,
            shift_encoding,
            unsafe_extract,
            RecodingScheme,
            Copyable,
            TwoToFour,
            FourToTwo,
            AsciiEncode,
            GenericRecoding
        """
        eval(Meta.parse(str))
    end
end

# Kmers.jl is tightly coupled to BioSequences and relies on much of its internals.
# Hence, we do not care about carefully importing specific symbols
using BioSequences
using BioSymbols: BioSymbol

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
const BitInteger =
    Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Int128, UInt128}

include("tuple_bitflipping.jl")
include("kmer.jl")
include("construction.jl")
include("indexing.jl")
include("transformations.jl")
include("revtrans.jl")
include("counting.jl")

include("iterators/common.jl")
include("iterators/FwKmers.jl")
include("iterators/CanonicalKmers.jl")
include("iterators/UnambiguousKmers.jl")
include("iterators/SpacedKmers.jl")

if !isdefined(Base, :get_extension)
    include("../ext/StringViewsExt.jl")
end

end # module
