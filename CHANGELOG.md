# Changelog
All notable changes to this project will be documented in this file.

## [1.1.0]
### Additions
* New exported function `as_integer` gets the integer representation of a kmer,
  if it's 128 bits or lower. Note that the precise integer is subject to change
  in minor versions of Kmers - see the docstring of the function for details.
* Added the ability to sample random kmers with `rand(::Type{<:Kmer})`
* Added the ability to sample symbols from a kmer with `rand(::Kmer)`
* Added optimised `count(isGC, ::Kmer{<:NucleicAcidAlphabet{2}})`.

## [1.0.0]
Complete overhaul of the API, including new types and new functions.
The changes are all-encompassing, so listing them out is pointless. For all
intents and purposes, Kmers v1.0.0 is a different package.
