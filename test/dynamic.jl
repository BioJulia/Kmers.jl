using Test
using Kmers
using BioSequences

@testset "Construction" begin
    for s in [
            dna"ATGTCGTTAGT",
            dna"",
            dna"YATGC-ATGwTCTDV",
            rna"AUGUCGAGUGUAUGC",
            rna"AUGCWSYGANN--CA",
            aa"KWOP",
            aa"",
            aa"TYAPC",

        ]
        for i in (s, string(s), (i for i in s))
            m = DynamicKmer{typeof(Alphabet(s)), UInt64}(i)
            @test length(m) == length(i)
            @test m == s
        end
    end

    # Two to four alphabet
    for s in Any[
        dna"ATGCTGTGACCA",
        dna"ATGTCGA",
        dna"",
    ]
        for A in Any[DNAAlphabet, RNAAlphabet]
            for (srcB, dstB) in [(2, 4), (4, 2)]
                src = LongSequence{A{srcB}}(s)
                dst = DynamicKmer{A{dstB}, UInt64}(src)

                @test src == dst
            end
        end
    end
end

@testset "Misc" begin
    s = dna"ATGCTGAC"
    m = DynamicDNAKmer{UInt32}(s)
    @test reverse_complement(m) == typeof(m)(reverse_complement(s))
end