using Test
using Kmers
using BioSequences

@testset "Construction" begin
    @testset "Same alphabet, and ASCII alphabet" begin
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
    end

    @testset "Two to four bit alphabet" begin
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

    @testset "Four to two bit alphabet" begin
        for s in [dna"TAGCTGAC", dna"ATGCTA", dna""]
            for A in [DNAAlphabet{2}, RNAAlphabet{2}]
                m = DynamicKmer{A, UInt64}(LongSequence{DNAAlphabet{4}}(s))
                @test m == LongSequence{A}(s)
            end
        end
    end

    @testset "Generic alphabet" begin
        for s in ["HE", "", "中Å!"]
            m = DynamicKmer{CharAlphabet, UInt128}(s)
            @test length(m) == length(s)
            @test string(m) == s
        end
    end

    @testset "From Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(s)
            dkmer = DynamicKmer{typeof(Alphabet(s)), UInt64}(kmer)
            @test dkmer == s
            @test length(dkmer) == length(kmer)
        end
    end

    @testset "To Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            dkmer = DynamicKmer{typeof(Alphabet(s)), UInt64}(s)
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(dkmer)
            @test length(kmer) == length(dkmer)
            @test string(dkmer) == string(kmer)
            @test_throws MethodError kmer == dkmer
        end
    end

    @testset "Capacity limits" begin
        # Test that exceeding capacity throws
        @test_throws Exception DynamicDNAKmer{UInt32}(dna"T"^30)
        @test_throws Exception DynamicAAKmer{UInt32}(aa"A"^8)
    end
end

@testset "Indexing and iteration" begin
    @testset "Scalar indexing" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGAC")
        @test m[1] == DNA_T
        @test m[3] == DNA_G
        @test m[8] == DNA_C
        @test first(m) == DNA_T
        @test last(m) == DNA_C

        @test_throws BoundsError m[0]
        @test_throws BoundsError m[9]
    end

    @testset "Range indexing" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGAC")
        @test m[1:3] == DynamicDNAKmer{UInt64}(dna"TAG")
        @test m[2:5] == DynamicDNAKmer{UInt64}(dna"AGCT")
        @test m[6:8] == DynamicDNAKmer{UInt64}(dna"GAC")
        @test m[1:0] == DynamicDNAKmer{UInt64}(dna"")

        @test_throws BoundsError m[0:3]
        @test_throws BoundsError m[6:9]
    end

    @testset "Iteration" begin
        s = dna"TAGCTGAC"
        m = DynamicDNAKmer{UInt64}(s)
        @test collect(m) == collect(s)
    end
end

@testset "Comparison and equality" begin
    @testset "Equality" begin
        @test DynamicDNAKmer{UInt64}(dna"TAG") == DynamicDNAKmer{UInt64}(dna"TAG")
        @test DynamicDNAKmer{UInt64}(dna"TAG") != DynamicDNAKmer{UInt64}(dna"TAC")
        @test DynamicDNAKmer{UInt64}(dna"") == DynamicDNAKmer{UInt64}(dna"")
    end

    @testset "Ordering" begin
        @test DynamicDNAKmer{UInt64}(dna"TAG") < DynamicDNAKmer{UInt64}(dna"TGA")
        @test DynamicDNAKmer{UInt64}(dna"AAA") < DynamicDNAKmer{UInt64}(dna"AAC")
        @test DynamicDNAKmer{UInt64}(dna"TAG") > DynamicDNAKmer{UInt64}(dna"TAC")
    end

    @testset "Comparison of different lengths" begin
        @test DynamicDNAKmer{UInt64}(dna"TAG") < DynamicDNAKmer{UInt64}(dna"TAGA")
        @test DynamicDNAKmer{UInt64}(dna"TAGA") > DynamicDNAKmer{UInt64}(dna"TAG")
    end
end

@testset "Integer conversion" begin
    @testset "as_integer" begin
        m1 = DynamicDNAKmer{UInt64}(dna"TAG")
        u1 = as_integer(m1)
        @test u1 isa Unsigned

        m2 = DynamicDNAKmer{UInt64}(dna"TAC")
        u2 = as_integer(m2)
        @test u1 != u2

        # Empty kmer
        @test as_integer(DynamicDNAKmer{UInt64}(dna"")) == 0
    end

    @testset "from_integer" begin
        for s in [dna"TAG", dna"TAGCTGA", dna"ATGCTAGC"]
            m = DynamicDNAKmer{UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m == m2
        end

        # Error on exceeding capacity
        @test_throws Exception from_integer(DynamicDNAKmer{UInt32}, UInt32(0), 30)
    end

    @testset "Round-trip conversion" begin
        for s in [dna"TAGCTGA", rna"UGCUGA", aa"PLKWM"]
            m = DynamicKmer{typeof(Alphabet(s)), UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m === m2
        end
    end
end

@testset "Hashing" begin
    m1 = DynamicDNAKmer{UInt64}(dna"TAG")
    m2 = DynamicDNAKmer{UInt64}(dna"TAG")
    m3 = DynamicDNAKmer{UInt64}(dna"TAC")

    # Same kmers hash to same value
    @test hash(m1) == hash(m2)

    # Different kmers likely hash to different values (not guaranteed but likely)
    @test hash(m1) != hash(m3)

    # Hash with seed
    h1 = hash(m1, UInt(123))
    h2 = hash(m1, UInt(456))
    @test h1 != h2

    # Length is part of hash
    m1 = DynamicDNAKmer{UInt64}(dna"TCA")
    m2 = DynamicDNAKmer{UInt64}(dna"TC")
    @test hash(m1) != hash(m2)
end

@testset "Biological operations" begin
    @testset "Reverse" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGA")
        @test reverse(m) == DynamicDNAKmer{UInt64}(dna"AGTCGAT")
        @test reverse(DynamicDNAKmer{UInt64}(dna"")) == DynamicDNAKmer{UInt64}(dna"")
    end

    @testset "Complement" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGA")
        @test complement(m) == DynamicDNAKmer{UInt64}(dna"ATCGACT")

        m2 = DynamicRNAKmer{UInt64}(rna"UAGCUGA")
        @test complement(m2) == DynamicRNAKmer{UInt64}(rna"AUCGACU")
    end

    @testset "Reverse complement" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGA")
        @test reverse_complement(m) == DynamicDNAKmer{UInt64}(dna"TCAGCTA")
    end

    @testset "Canonical" begin
        m1 = DynamicDNAKmer{UInt64}(dna"TAGCTGA")
        m2 = DynamicDNAKmer{UInt64}(dna"TCAGCTA")
        @test canonical(m1) == canonical(m2)
        @test iscanonical(DynamicDNAKmer{UInt64}(dna"AATT"))
        @test iscanonical(DynamicDNAKmer{UInt64}(dna"TTAA"))
        @test iscanonical(empty(m1))
        @test !iscanonical(DynamicDNAKmer{UInt64}(dna"TGGA"))
    end
end

@testset "Counting" begin
    @testset "Count GC" begin
        @test count(isGC, DynamicDNAKmer{UInt64}(dna"TATCGGAGA")) == 4
        @test count(isGC, DynamicDNAKmer{UInt64}(dna"TATATATAAAAA")) == 0
        @test count(isGC, DynamicDNAKmer{UInt64}(dna"")) == 0

        @test count(isGC, DynamicRNAKmer{UInt64}(rna"AUGUCGUAG")) == 4
    end

    @testset "Count symbols" begin
        m = DynamicDNAKmer{UInt64}(dna"TAGCTGA")
        @test count(==(DNA_A), m) == 2
        @test count(==(DNA_T), m) == 2
        @test count(==(DNA_G), m) == 2
        @test count(==(DNA_C), m) == 1
    end
end

@testset "shift_encoding" begin
    m = DynamicDNAKmer{UInt32}(dna"TAGA")
    enc = UInt32(BioSequences.encode(DNAAlphabet{2}(), DNA_C))
    m2 = Kmers.shift_encoding(m, enc)
    @test m2 == DynamicDNAKmer{UInt32}(dna"AGAC")
    @test length(m2) == length(m)
end

@testset "Mixed integer types" begin
    for U in [UInt8, UInt16, UInt32, UInt64, UInt128]
        s = dna"TAG"
        m = DynamicKmer{DNAAlphabet{2}, U}(s)
        @test m == s
        @test length(m) == 3
    end
end
