# <img src="./sticker.svg" width="30%" align="right" /> Kmers

[![Latest Release](https://img.shields.io/github/release/BioJulia/Kmers.jl.svg)](https://github.com/BioJulia/Kmers.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/Kmers.jl/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/Kmers.jl/stable)
[![Pkg Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


## Description

Kmers provides a specialised concrete `BioSequence` subtype, optimised for
representing short immutable sequences called kmers: contiguous sub-strings of k
nucleotides of some reference sequence.

They are used extensively in bioinformatic analyses as an informational unit.
This concept was popularised by short read assemblers. 
Analyses within the kmer space benefit from a simple formulation of the sampling
problem and direct in-hash comparisons.

Kmers provides the type representing kmers as well as the implementations of
the APIs specified by the
[`BioSequences.jl`](https://github.com/BioJulia/BioSequences.jl) package.

## Installation

You can install BioSequences from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```julia
add Kmers
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

Kmers is tested against Julia `1.X` on Linux, OS X, and Windows.

[![Unit tests](https://github.com/BioJulia/Kmers.jl/workflows/Unit%20tests/badge.svg?branch=master)](https://github.com/BioJulia/Kmers.jl/actions?query=workflow%3A%22Unit+tests%22+branch%3Amaster)
[![Documentation](https://github.com/BioJulia/Kmers.jl/workflows/Documentation/badge.svg?branch=master)](https://github.com/BioJulia/BioKmers.jl/actions?query=workflow%3ADocumentation+branch%3Amaster)
[![](https://codecov.io/gh/BioJulia/Kmers.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BioJulia/Kmers.jl)


## Contributing

We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.


## Questions?

If you have a question about contributing or using BioJulia software, come
on over and chat to us on [Gitter](https://gitter.im/BioJulia/General), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).
