```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
end
```
## Kmer composition
In metagenomics, sequences are often summarized by counting the occurrence of
all k-mers of a given length in a sequence.
For example, for K=4, there are 4^4 = 256 possible DNA 4-mers.
If these counts are ordered, the composition can be represented by a length 256
vector.

Vector similarity operations (e.g. cosine distance) can then be used as an
approximate proxy for phylogenetic distance.

In the example below, we exploit that:
* A `DNAKmer{4}`'s data is a single-element tuple, which
  stores the sequence in the 8 lower bits.
* The `encoded_data` function will return this tuple.

```jldoctest; output=false
using BioSequences, FASTX, Kmers
using BioSequences: encoded_data

function composition(record::FASTARecord)
    counts = zeros(UInt32, 256)
    frequencies = zeros(Float32, 256)
    for kmer in FwDNAMers{4}(sequence(record))
        @inbounds counts[only(encoded_data(kmer)) + 1] += 1
    end
    factor = 1 / sum(counts; init=zero(eltype(counts)))
    for i in eachindex(counts, frequencies)
        frequencies[i] = counts[i] * factor
    end
    frequencies
end

# Make two FASTA records - could be from an assembly
recs = [FASTARecord(string(i), randdnaseq(10000)) for i in "AB"]

# Compute the 2-norm difference and verify it's in [0, 2].
(comp_a, comp_b) = map(composition, recs)
comp_distance = sum((comp_a .- comp_b).^2)
println(0.0 ≤ comp_distance ≤ 2.0)

# output
true

```