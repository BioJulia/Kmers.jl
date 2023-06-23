@testset "Comparisons" begin
    @testset "Equality" begin
        function check_seq_kmer_equality(len)
            a = DNAKmer(random_dna_kmer(len))
            b = LongDNA{4}(a)
            c = LongDNA{2}(a)
            return a == b == c && c == b == a
        end

        for len in [1, 10, 32, 64, 128]
            @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
        end

        # True negatives
        @test DNAKmer("ACG") != RNAKmer("ACG")
        @test DNAKmer("T") != RNAKmer("U")
        @test DNAKmer("AC") != DNAKmer("AG")
        @test RNAKmer("AC") != RNAKmer("AG")
        @test AAKmer("MV") != AAKmer("NM")

        @test DNAKmer("ACG") != rna"ACG"
        @test DNAKmer("T") != rna"U"
        @test DNAKmer("AC") != dna"AG"
        @test RNAKmer("AC") != rna"AG"
        @test AAKmer("MV") != aa"NM"

        @test rna"ACG" != DNAKmer("ACG")
        @test rna"U" != DNAKmer("T")
        @test dna"AG" != DNAKmer("AC")
        @test rna"AG" != RNAKmer("AC")
        @test aa"MV" != AAKmer("NM")
    end

    @testset "Inequality" begin
        for len in [1, 10, 32, 64]
            if len <= 32
                @test isless(DNAKmer{1}((UInt64(0),)), DNAKmer{1}((UInt64(1),)))
                @test !isless(DNAKmer{1}((UInt64(0),)), DNAKmer{1}((UInt64(0),)))
                @test !isless(DNAKmer{1}((UInt64(1),)), DNAKmer{1}((UInt64(0),)))

                @test isless(RNAKmer{1}((UInt64(0),)), RNAKmer{1}((UInt64(1),)))
                @test !isless(RNAKmer{1}((UInt64(0),)), RNAKmer{1}((UInt64(0),)))
                @test !isless(RNAKmer{1}((UInt64(1),)), RNAKmer{1}((UInt64(0),)))
            end
        end
    end

    @testset "Hash" begin
        kmers = map(DNAKmer, ["AAAA", "AACT", "ACGT", "TGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end

        kmers = map(RNAKmer, ["AAAA", "AACU", "ACGU", "UGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end

        kmers = map(AAKmer, ["AMVK", "FPST", "QEGH", "ARND"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
    end
end
