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

    @testset "To Kmer with N=2 (>64 coding bits)" begin
        # For 2-bit DNA: need >32 bases for >64 coding bits
        s_dna = dna"TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"  # 33 bases = 66 bits
        dkmer_dna = DynamicKmer{DNAAlphabet{2}, UInt128}(s_dna)
        kmer_dna = Kmer{DNAAlphabet{2}, 33, 2}(dkmer_dna)
        @test length(kmer_dna) == length(dkmer_dna)
        @test string(dkmer_dna) == string(kmer_dna)

        # For 2-bit RNA
        s_rna = rna"UGCUGAUGCUGAUGCUGAUGCUGAUGCUGAUGA"  # 33 bases = 66 bits
        dkmer_rna = DynamicKmer{RNAAlphabet{2}, UInt128}(s_rna)
        kmer_rna = Kmer{RNAAlphabet{2}, 33, 2}(dkmer_rna)
        @test length(kmer_rna) == length(dkmer_rna)
        @test string(dkmer_rna) == string(kmer_rna)

        # For 8-bit amino acids: need >8 bases for >64 coding bits
        s_aa = aa"KWOPPLKWM"  # 9 bases = 72 bits
        dkmer_aa = DynamicKmer{AminoAcidAlphabet, UInt128}(s_aa)
        kmer_aa = Kmer{AminoAcidAlphabet, 9, 2}(dkmer_aa)
        @test length(kmer_aa) == length(dkmer_aa)
        @test string(dkmer_aa) == string(kmer_aa)

        # Test error on length mismatch
        dkmer = DynamicDNAKmer{UInt64}(dna"TAG")
        @test_throws Exception Kmer{DNAAlphabet{2}, 5, 1}(dkmer)
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

    @testset "cmp function" begin
        # Test cmp returns -1, 0, or 1
        @test cmp(DynamicDNAKmer{UInt64}(dna"TAG"), DynamicDNAKmer{UInt64}(dna"TAG")) == 0
        @test cmp(DynamicDNAKmer{UInt64}(dna"TAG"), DynamicDNAKmer{UInt64}(dna"TGA")) < 0
        @test cmp(DynamicDNAKmer{UInt64}(dna"TGA"), DynamicDNAKmer{UInt64}(dna"TAG")) > 0

        # Test with RNA
        @test cmp(DynamicRNAKmer{UInt64}(rna"UAG"), DynamicRNAKmer{UInt64}(rna"UAG")) == 0
        @test cmp(DynamicRNAKmer{UInt64}(rna"UAG"), DynamicRNAKmer{UInt64}(rna"UGA")) < 0

        # Test with different lengths
        @test cmp(DynamicDNAKmer{UInt64}(dna"TAG"), DynamicDNAKmer{UInt64}(dna"TAGA")) < 0
        @test cmp(DynamicDNAKmer{UInt64}(dna"TAGA"), DynamicDNAKmer{UInt64}(dna"TAG")) > 0
    end

    @testset "Comparison across different integer types" begin
        s1 = dna"TAG"
        s2 = dna"TAC"
        s3 = dna"TAGA"

        m32_1 = DynamicDNAKmer{UInt32}(s1)
        m64_1 = DynamicDNAKmer{UInt64}(s1)
        m32_2 = DynamicDNAKmer{UInt32}(s2)
        m64_2 = DynamicDNAKmer{UInt64}(s2)

        # Test equality across types
        @test m32_1 == m64_1
        @test m32_2 == m64_2
        @test m32_1 != m64_2

        # Test ordering across types
        @test m32_2 < m64_1  # TAC < TAG
        @test m64_1 > m32_2  # TAG > TAC

        # Test cmp across types
        @test cmp(m32_1, m64_1) == 0
        @test cmp(m32_2, m64_1) < 0
        @test cmp(m64_1, m32_2) > 0

        # Test with different lengths
        m32_3 = DynamicDNAKmer{UInt32}(s3)
        @test m32_1 < m32_3
        @test m64_1 < m32_3
        @test cmp(m32_1, m32_3) < 0

        # Also for AA types
        s1 = aa"KPRCRLF"
        s2 = aa"KPRCRLFAAA"

        m64_1 = DynamicAAKmer{UInt64}(s1)
        m128_1 = DynamicAAKmer{UInt128}(s1)
        m128_2 = DynamicAAKmer{UInt128}(s2)

        @test m64_1 == m128_1
        @test cmp(m64_1, m128_1) == 0
        @test m128_1 != m128_2
        @test m128_1 < m128_2
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

    # Length must be part of hash: TCA and TC have identical coding bits
    # (since A encodes to 00, which looks like padding), but different lengths.
    # They must hash differently despite having the same as_integer representation.
    m1 = DynamicDNAKmer{UInt64}(dna"TCA")
    m2 = DynamicDNAKmer{UInt64}(dna"TC")
    @test hash(m1) != hash(m2)

    @testset "fx_hash" begin
        m1 = DynamicDNAKmer{UInt64}(dna"TAG")
        m2 = DynamicDNAKmer{UInt64}(dna"TAG")
        m3 = DynamicDNAKmer{UInt64}(dna"TAC")

        # Same kmers should produce same fx_hash
        @test Kmers.fx_hash(m1, UInt64(0)) == Kmers.fx_hash(m2, UInt64(0))

        # Different kmers should produce different fx_hash values
        @test Kmers.fx_hash(m1, UInt64(0)) != Kmers.fx_hash(m3, UInt64(0))

        # Different seeds should produce different hashes
        @test Kmers.fx_hash(m1, UInt64(123)) != Kmers.fx_hash(m1, UInt64(456))

        # Length must be part of hash: TCA and TC have identical coding bits
        # (since A encodes to 00, which looks like padding), but different lengths.
        # They must hash differently despite having the same as_integer representation.
        m1 = DynamicDNAKmer{UInt64}(dna"TCA")
        m2 = DynamicDNAKmer{UInt64}(dna"TC")
        @test Kmers.fx_hash(m1, UInt64(0)) != Kmers.fx_hash(m2, UInt64(0))
    end
end

@testset "Biological operations" begin
    @testset "Reverse" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = DynamicDNAKmer{UInt64}(s)
        @test reverse(m) == DynamicDNAKmer{UInt64}(reverse(s))

        # Test empty sequence
        @test reverse(DynamicDNAKmer{UInt64}(dna"")) == DynamicDNAKmer{UInt64}(dna"")

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = DynamicRNAKmer{UInt64}(s_rna)
        @test reverse(m_rna) == DynamicRNAKmer{UInt64}(reverse(s_rna))

        # Test with 4-bit DNA
        s_4bit = dna"TAGCTGA"
        m_4bit = DynamicKmer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse(m_4bit) == DynamicKmer{DNAAlphabet{4}, UInt64}(reverse(s_4bit))

        # Test with amino acids
        s_aa = aa"KWOP"
        m_aa = DynamicAAKmer{UInt64}(s_aa)
        @test reverse(m_aa) == DynamicAAKmer{UInt64}(reverse(s_aa))
    end

    @testset "Complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = DynamicDNAKmer{UInt64}(s)
        @test complement(m) == DynamicDNAKmer{UInt64}(complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = DynamicRNAKmer{UInt64}(s_rna)
        @test complement(m_rna) == DynamicRNAKmer{UInt64}(complement(s_rna))

        # Test with 4-bit DNA (includes ambiguous bases)
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGAW")  # W = A or T
        m_4bit = DynamicKmer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test complement(m_4bit) == DynamicKmer{DNAAlphabet{4}, UInt64}(complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGAW")  # W = A or U
        m_rna_4bit = DynamicKmer{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test complement(m_rna_4bit) == DynamicKmer{RNAAlphabet{4}, UInt64}(complement(s_rna_4bit))
    end

    @testset "Reverse complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = DynamicDNAKmer{UInt64}(s)
        @test reverse_complement(m) == DynamicDNAKmer{UInt64}(reverse_complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = DynamicRNAKmer{UInt64}(s_rna)
        @test reverse_complement(m_rna) == DynamicRNAKmer{UInt64}(reverse_complement(s_rna))

        # Test with 4-bit DNA
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGA")
        m_4bit = DynamicKmer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse_complement(m_4bit) == DynamicKmer{DNAAlphabet{4}, UInt64}(reverse_complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGA")
        m_rna_4bit = DynamicKmer{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test reverse_complement(m_rna_4bit) == DynamicKmer{RNAAlphabet{4}, UInt64}(reverse_complement(s_rna_4bit))
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

        m = DynamicAAKmer{UInt128}(aa"WLAKWVMARQKW")
        @test count(==(AA_W), m) == 3
        @test count(==(AA_Q), m) == 1
        @test count(==(AA_A), m) == 2
        @test count(==(AA_C), m) == 0
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

@testset "Misc" begin
    d = DynamicAAKmer{UInt32}("WPK")
    @test only([d]') === d
end
