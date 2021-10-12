# Kmer types

Bioinformatic analyses make extensive use of kmers.
Kmers are contiguous sub-strings of k nucleotides of some reference sequence. 

They are used extensively in bioinformatic analyses as an informational unit.
This concept popularised by short read assemblers. 
Analyses within the kmer space benefit from a simple formulation of the sampling
problem and direct in-hash comparisons.

BioSequences provides the following types to represent Kmers.

```@docs
Kmer
```

The following aliases are also defined:

```@docs
DNAKmer
DNAKmer27
DNAKmer31
DNAKmer63
RNAKmer
RNAKmer27
RNAKmer31
RNAKmer63
```


### Skipmers

For some analyses, the contiguous nature of kmers imposes limitations.
A single base difference, due to real biological variation or a sequencing error,
affects all k-mers crossing that position thus impeding direct analyses by identity.
Also, given the strong interdependence of local sequence, contiguous sections
capture less information about genome structure, and so they are more affected by
sequence repetition. 

Skipmers are a generalisation of the concept of a kmer.
They are created using a cyclic pattern of used-and-skipped positions which
achieves increased entropy and tolerance to nucleotide substitution differences
by following some simple rules.

Skipmers preserve many of the elegant properties of kmers such as reverse
complementability and existence of a canonical representation.
Also, using cycles of three greatly increases the power of direct intersection
between the genomes of different organisms by grouping together the more conserved 
nucleotides of protein-coding regions.

BioSequences currently does not provide a separate type for skipmers, they are
represented using `Mer` and `BigMer` as their representation as a short immutable
sequence encoded in an unsigned integer is the same.
The distinction lies in how they are generated.

#### Skipmer generation

A skipmer is a simple cyclic q-gram that includes _m_ out of every _n_ bases
until a total of _k_ bases is reached. 

This is illustrated in the figure below (from
[this paper](https://www.biorxiv.org/content/biorxiv/early/2017/08/23/179960.full.pdf).):

![skipmer-fig](skipmers.png)

To maintain cyclic properties and the existence of the reverse-complement as a
skipmer defined by the same function, _k_ should be a multiple of _m_.

This also enables the existence of a canonical representation for each skipmer,
defined as the lexicographically smaller of the forward and reverse-complement 
representations.

Defining _m_, _n_ and _k_ fixes a value for _S_, the total span of the skipmer,
given by: 

```math
S = n * (\frac{k}{m} - 1) + m
```

To see how to iterate over skipmers cf. kmers, see the Iteration section
of the manual.

