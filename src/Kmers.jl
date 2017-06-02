# Kmers.jl
# ========
#
# Module for kmers and short sequence fragments.
#
# This file is a part of the Kmers.jl, a package in the BioJulia ecosystem.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

__precompile__()

module Kmers

import BioSymbols: NucleicAcid, DNA, RNA

include("kmer.jl")

end # module
