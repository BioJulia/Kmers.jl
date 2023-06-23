```@meta
CurrentModule = Kmers
DocTestSetup = quote
    using BioSequences
    using Test
    using Kmers
    using FASTX
    using MinHash
end
```
## MinHash
The MinHash algorithm is used in tools such as
[Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)
and [sourmash](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/)
to quickly compute approximate similarities of genomes, collections of genomes, or collections of reads.

```jldoctest; filter = r"^\d+ MB/s$" => s"***"
using BioSequences, MinHash, FASTX, Kmers

# Write 25 sequences of length 20 to a buffer.
# Try changing this to length 4 million!
buffer = IOBuffer()
writer = FASTAWriter(buffer)
n_bytes = sum(1:25) do genome
    rec = FASTARecord("seq_$(genome)", randdnaseq(20))
    write(writer, rec)
end
flush(writer)

# Time minhashing the 50 genomes
timing = @timed FASTAReader(seekstart(buffer); copy=false) do reader
    map(reader) do record
        seq = codeunits(sequence(record))
        sketch(fx_hash, CanonicalDNAMers{16}(sequence(record)), 1000)
    end
end
println(round(Int, n_bytes / (timing.time * 1e6)), " MB/s")

# output

200 MB/s
```
