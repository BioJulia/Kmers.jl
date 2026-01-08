using Test
using Kmers
using BioSequences
using BitIntegers

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

        # Test macro equivalence with constructor
        @test dmer"ATGTCGTTAGT"d == DynamicDNAKmer{UInt64}(dna"ATGTCGTTAGT")
        @test dmer"AUGUCGAGUGUAUGC"r == DynamicRNAKmer{UInt64}(rna"AUGUCGAGUGUAUGC")
        @test dmer"KWOP"a == DynamicAAKmer{UInt64}(aa"KWOP")

        d = DynamicDNAKmer{UInt}(dna"ATGTCGTTAGT")
        @test DynamicRNAKmer(d) == d
        @test DynamicRNAKmer{UInt64}(d) == d

        # From a large kmer
        m = mer"TAGTGCTGTAGTAGTGCTGTATGATGTCTGCATGC"d
        dm = DynamicDNAKmer{UInt256}(m)
        @test LongSequence(m) == LongSequence(dm)

        # Using large bitintegers - 240 2-bit nt = 480 bits should fit in
        # 512 bits with plenty room for the runtime length
        seq = randdnaseq(240)
        dm = DynamicDNAKmer{UInt512}(seq)
        @test string(seq) == string(dm)
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
            m = DynamicKmer{CharAlphabet, UInt512}(s)
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

    @testset "From DynamicKmer" begin
        # Same backing type, same alphabet
        s = dna"TAGCTGA"
        d1 = dmer"TAGCTGA"d
        d2 = DynamicDNAKmer{UInt64}(d1)
        @test d2 == d1
        @test d2 === d1  # Should be identical
        @test length(d2) == length(d1)

        # Different backing type, same alphabet - widening
        d32 = DynamicDNAKmer{UInt32}(dna"TAGC")
        d64 = DynamicDNAKmer{UInt64}(d32)
        @test d64 == d32
        @test length(d64) == length(d32)
        @test d64 == dna"TAGC"

        # Different backing type, same alphabet - narrowing
        d128 = DynamicDNAKmer{UInt128}(dna"TAGC")
        d32_narrow = DynamicDNAKmer{UInt32}(d128)
        @test d32_narrow == d128
        @test length(d32_narrow) == length(d128)
        @test d32_narrow == dna"TAGC"

        # Same alphabet family (DNA/RNA), different alphabets
        d_dna = dmer"ATGTCGTTAGT"d
        d_rna = DynamicRNAKmer{UInt64}(d_dna)
        @test d_rna == d_dna
        @test length(d_rna) == length(d_dna)

        # Different backing type AND different alphabet
        d_dna32 = DynamicDNAKmer{UInt32}(dna"TAGC")
        d_rna64 = DynamicRNAKmer{UInt64}(d_dna32)
        @test d_rna64 == d_dna32
        @test length(d_rna64) == length(d_dna32)

        # Test with amino acids
        aa_d64 = dmer"KWOP"a
        aa_d128 = DynamicAAKmer{UInt128}(aa_d64)
        @test aa_d128 == aa_d64
        @test length(aa_d128) == length(aa_d64)
        aa_d256 = DynamicAAKmer{UInt256}(aa_d64)
        @test aa_d256 == aa_d64
        @test length(aa_d256) == length(aa_d64)

        # Test with empty kmers
        empty_d32 = DynamicDNAKmer{UInt32}(dna"")
        empty_d64 = DynamicDNAKmer{UInt64}(empty_d32)
        @test empty_d64 == empty_d32
        @test length(empty_d64) == 0

        # Test with 4-bit alphabets
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCN")
        d_4bit_32 = DynamicKmer{DNAAlphabet{4}, UInt32}(s_4bit)
        d_4bit_64 = DynamicKmer{DNAAlphabet{4}, UInt64}(d_4bit_32)
        @test d_4bit_64 == d_4bit_32
        @test length(d_4bit_64) == length(d_4bit_32)
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
        dkmer_dna = DynamicKmer{DNAAlphabet{2}, UInt256}(s_dna)
        kmer_dna = Kmer{DNAAlphabet{2}, 33}(dkmer_dna)
        @test length(kmer_dna) == length(dkmer_dna)
        @test string(dkmer_dna) == string(kmer_dna)

        # For 2-bit RNA
        s_rna = rna"UGCUGAUGCUGAUGCUGAUGCUGAUGCUGAUGA"  # 33 bases = 66 bits
        dkmer_rna = DynamicKmer{RNAAlphabet{2}, UInt256}(s_rna)
        kmer_rna = Kmer{RNAAlphabet{2}, 33}(dkmer_rna)
        @test length(kmer_rna) == length(dkmer_rna)
        @test string(dkmer_rna) == string(kmer_rna)

        # For 8-bit amino acids: need >8 bases for >64 coding bits
        s_aa = aa"KWOPPLKWM"  # 9 bases = 72 bits
        dkmer_aa = DynamicKmer{AminoAcidAlphabet, UInt256}(s_aa)
        kmer_aa = Kmer{AminoAcidAlphabet, 9}(dkmer_aa)
        @test length(kmer_aa) == length(dkmer_aa)
        @test string(dkmer_aa) == string(kmer_aa)

        # Test error on length mismatch
        dkmer = dmer"TAG"d
        @test_throws Exception Kmer{DNAAlphabet{2}, 5}(dkmer)
    end

    @testset "Capacity limits" begin
        # Test that exceeding capacity throws
        @test_throws Exception DynamicDNAKmer{UInt32}(dna"T"^30)
        @test_throws Exception DynamicAAKmer{UInt32}(aa"A"^8)
    end
end

@testset "Indexing and iteration" begin
    @testset "Scalar indexing" begin
        m = dmer"TAGCTGAC"d
        @test m[1] == DNA_T
        @test m[3] == DNA_G
        @test m[8] == DNA_C
        @test first(m) == DNA_T
        @test last(m) == DNA_C

        @test_throws BoundsError m[0]
        @test_throws BoundsError m[9]
    end

    @testset "Range indexing" begin
        m = dmer"TAGCTGAC"d
        @test m[1:3] == dmer"TAG"d
        @test m[2:5] == dmer"AGCT"d
        @test m[6:8] == dmer"GAC"d
        @test m[1:0] == dmer""d

        @test_throws BoundsError m[0:3]
        @test_throws BoundsError m[6:9]
    end

    @testset "Iteration" begin
        s = dna"TAGCTGAC"
        m = dmer"TAGCTGAC"d
        @test collect(m) == collect(s)
    end
end

@testset "Comparison and equality" begin
    @testset "Equality" begin
        @test dmer"TAG"d == dmer"TAG"d
        @test dmer"TAG"d != dmer"TAC"d
        @test dmer""d == dmer""d
    end

    @testset "Ordering" begin
        @test dmer"TAG"d < dmer"TGA"d
        @test dmer"AAA"d < dmer"AAC"d
        @test dmer"TAG"d > dmer"TAC"d
    end

    @testset "Comparison of different lengths" begin
        @test dmer"TAG"d < dmer"TAGA"d
        @test dmer"TAGA"d > dmer"TAG"d
    end

    @testset "cmp function" begin
        # Test cmp returns -1, 0, or 1
        @test cmp(dmer"TAG"d, dmer"TAG"d) == 0
        @test cmp(dmer"TAG"d, dmer"TGA"d) < 0
        @test cmp(dmer"TGA"d, dmer"TAG"d) > 0

        # Test with RNA
        @test cmp(dmer"UAG"r, dmer"UAG"r) == 0
        @test cmp(dmer"UAG"r, dmer"UGA"r) < 0

        # Test with different lengths
        @test cmp(dmer"TAG"d, dmer"TAGA"d) < 0
        @test cmp(dmer"TAGA"d, dmer"TAG"d) > 0
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
        m256_1 = DynamicAAKmer{UInt256}(s1)

        @test m64_1 == m128_1
        @test cmp(m64_1, m128_1) == 0
        @test m128_1 != m128_2
        @test m128_1 < m128_2
        @test m64_1 == m256_1
        @test cmp(m256_1, m128_1) == 0
    end
end

@testset "Integer conversion" begin
    @testset "as_integer" begin
        m1 = dmer"TAG"d
        u1 = as_integer(m1)
        @test u1 isa Unsigned

        m2 = dmer"TAC"d
        u2 = as_integer(m2)
        @test u1 != u2

        # Empty kmer
        @test as_integer(dmer""d) == 0
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
    m1 = dmer"TAG"d
    m2 = dmer"TAG"d
    m3 = dmer"TAC"d

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
    m1 = dmer"TCA"d
    m2 = dmer"TC"d
    @test hash(m1) != hash(m2)

    @testset "fx_hash" begin
        m1 = dmer"TAG"d
        m2 = dmer"TAG"d
        m3 = dmer"TAC"d

        # Same kmers should produce same fx_hash
        @test Kmers.fx_hash(m1, UInt64(0)) == Kmers.fx_hash(m2, UInt64(0))

        # Different kmers should produce different fx_hash values
        @test Kmers.fx_hash(m1, UInt64(0)) != Kmers.fx_hash(m3, UInt64(0))

        # Different seeds should produce different hashes
        @test Kmers.fx_hash(m1, UInt64(123)) != Kmers.fx_hash(m1, UInt64(456))

        # Length must be part of hash: TCA and TC have identical coding bits
        # (since A encodes to 00, which looks like padding), but different lengths.
        # They must hash differently despite having the same as_integer representation.
        m1 = dmer"TCA"d
        m2 = dmer"TC"d
        @test Kmers.fx_hash(m1, UInt64(0)) != Kmers.fx_hash(m2, UInt64(0))
    end
end

@testset "Biological operations" begin
    @testset "Reverse" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test reverse(m) == DynamicDNAKmer{UInt64}(reverse(s))

        # Test empty sequence
        @test reverse(dmer""d) == dmer""d

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse(m_rna) == DynamicRNAKmer{UInt64}(reverse(s_rna))

        # Test with 4-bit DNA
        s_4bit = dna"TAGCTGA"
        m_4bit = DynamicKmer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse(m_4bit) == DynamicKmer{DNAAlphabet{4}, UInt64}(reverse(s_4bit))

        # Test with amino acids
        s_aa = aa"KWOP"
        m_aa = dmer"KWOP"a
        @test reverse(m_aa) == DynamicAAKmer{UInt64}(reverse(s_aa))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGC"
        m_large = DynamicDNAKmer{UInt256}(s_large)
        @test reverse(m_large) == DynamicDNAKmer{UInt256}(reverse(s_large))
    end

    @testset "Complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test complement(m) == DynamicDNAKmer{UInt64}(complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
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
        m = dmer"TAGCTGA"d
        @test reverse_complement(m) == DynamicDNAKmer{UInt64}(reverse_complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse_complement(m_rna) == DynamicRNAKmer{UInt64}(reverse_complement(s_rna))

        # Test with 4-bit DNA
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGA")
        m_4bit = DynamicKmer{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse_complement(m_4bit) == DynamicKmer{DNAAlphabet{4}, UInt64}(reverse_complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGA")
        m_rna_4bit = DynamicKmer{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test reverse_complement(m_rna_4bit) == DynamicKmer{RNAAlphabet{4}, UInt64}(reverse_complement(s_rna_4bit))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGCTAGC"
        m_large = DynamicDNAKmer{UInt512}(s_large)
        @test reverse_complement(m_large) == DynamicDNAKmer{UInt512}(reverse_complement(s_large))
    end

    @testset "Canonical" begin
        m1 = dmer"TAGCTGA"d
        m2 = dmer"TCAGCTA"d
        @test canonical(m1) == canonical(m2) == m1
        @test iscanonical(dmer"AATT"d)
        @test iscanonical(dmer"TTAA"d)
        @test iscanonical(empty(m1))
        @test !iscanonical(dmer"TGGA"d)

        # Test with larger bit integers
        m_large = DynamicDNAKmer{UInt256}(dna"TAGCTAGCTAGC")
        m_rc = reverse_complement(m_large)
        @test canonical(m_large) == min(m_large, m_rc)
    end
end

@testset "Counting" begin
    @testset "Count GC" begin
        @test count(isGC, dmer"TATCGGAGA"d) == 4
        @test count(isGC, dmer"TATATATAAAAA"d) == 0
        @test count(isGC, dmer""d) == 0

        @test count(isGC, dmer"AUGUCGUAG"r) == 4
    end

    @testset "Count symbols" begin
        m = dmer"TAGCTGA"d
        @test count(==(DNA_A), m) == 2
        @test count(==(DNA_T), m) == 2
        @test count(==(DNA_G), m) == 2
        @test count(==(DNA_C), m) == 1

        m = DynamicAAKmer{UInt256}(aa"WLAKWVMARQKW")
        @test count(==(AA_W), m) == 3
        @test count(==(AA_Q), m) == 1
        @test count(==(AA_A), m) == 2
        @test count(==(AA_C), m) == 0

        # Test with even larger integers
        m_large = DynamicDNAKmer{UInt512}(dna"TAGCTAGCTAGCTAGC")
        @test count(==(DNA_T), m_large) == 4
        @test count(==(DNA_A), m_large) == 4
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
    for U in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512]
        s = dna"TAG"
        m = DynamicKmer{DNAAlphabet{2}, U}(s)
        @test m == s
        @test length(m) == 3
    end
end

@testset "push and push_first" begin
    for (dkmer_fn, longseq_fn) in Any[[push, push!], [push_first, pushfirst!]]
        for (dkmer, symbol) in Any[
                (dmer"TAGC"d, DNA_A),
                (dmer"AUGC"r, RNA_U),
                (dmer""d, DNA_T),
                (dmer""r, RNA_G),
                (DynamicDNAKmer{UInt32}(dna"TAG"), DNA_G),
                (DynamicRNAKmer{UInt32}(rna"AUG"), RNA_C),
                (DynamicAAKmer{UInt64}(aa"KWOP"), AA_L),
                (DynamicAAKmer{UInt256}(aa"MWP"), AA_K),
            ]
            # Apply operation to DynamicKmer
            result_dkmer = dkmer_fn(dkmer, symbol)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq, symbol)

            # Test equality
            @test result_dkmer == longseq
            @test length(result_dkmer) == length(longseq)
            @test length(result_dkmer) == length(dkmer) + 1

            # Test that result has correct type
            @test result_dkmer isa typeof(dkmer)

            # Test that original is unchanged (immutability)
            @test length(dkmer) == length(result_dkmer) - 1
        end
    end

    # Test type conversion (pushing DNA to RNA kmer and vice versa)
    @test push(dmer"ATGC"d, RNA_U) == dmer"ATGCT"d
    @test push(dmer"AUGC"r, DNA_T) == dmer"AUGCU"r
    @test push_first(dmer"ATGC"d, RNA_U) == dmer"TATGC"d
    @test push_first(dmer"AUGC"r, DNA_T) == dmer"UAUGC"r

    # Test character conversion
    @test push(dmer"TAG"d, 'C') == dmer"TAGC"d
    @test push_first(dmer"TAG"d, 'A') == dmer"ATAG"d

    # Test chaining operations
    @test push(push(dmer"TAG"d, DNA_C), DNA_A) == dmer"TAGCA"d
    @test push_first(push_first(dmer"TAG"d, DNA_C), DNA_A) == dmer"ACTAG"d

    # Test error when kmer is at max capacity
    # For UInt32 with 2-bit DNA: capacity = 14
    d_full = DynamicDNAKmer{UInt32}(dna"T"^14)
    @test_throws BoundsError push(d_full, DNA_A)
    @test_throws BoundsError push_first(d_full, DNA_A)

    # For UInt64 with 2-bit DNA: capacity = 29
    d64_full = DynamicDNAKmer{UInt64}(dna"A"^29)
    @test_throws BoundsError push(d64_full, DNA_T)
    @test_throws BoundsError push_first(d64_full, DNA_G)

    # For UInt32 with 8-bit AA: capacity = 3
    aa32_full = DynamicAAKmer{UInt32}(aa"WPK")
    @test_throws BoundsError push(aa32_full, AA_L)
    @test_throws BoundsError push_first(aa32_full, AA_M)

    # Verify we can push to capacity-1 without error
    d_almost_full = DynamicDNAKmer{UInt64}(dna"T"^28)
    @test length(push(d_almost_full, DNA_A)) == 29
    @test length(push_first(d_almost_full, DNA_G)) == 29

    # Test with larger bit integers
    d256 = DynamicDNAKmer{UInt256}(dna"TAGC")
    @test push(d256, DNA_G) == DynamicDNAKmer{UInt256}(dna"TAGCG")
    @test push_first(d256, DNA_A) == DynamicDNAKmer{UInt256}(dna"ATAGC")
end

@testset "pop and pop_first" begin
    for (dkmer_fn, longseq_fn) in Any[[pop, pop!], [pop_first, popfirst!]]
        for dkmer in Any[
                dmer"TAGCA"d,
                dmer"AUGCU"r,
                dmer"T"d,
                dmer"A"r,
                DynamicDNAKmer{UInt32}(dna"TAGG"),
                DynamicRNAKmer{UInt32}(rna"AUGC"),
                DynamicAAKmer{UInt64}(aa"KWOPL"),
                DynamicAAKmer{UInt128}(aa"MWPK"),
            ]
            # Apply operation to DynamicKmer
            result_dkmer = dkmer_fn(dkmer)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq)

            # Test equality
            @test result_dkmer == longseq
            @test length(result_dkmer) == length(longseq)
            @test length(result_dkmer) == length(dkmer) - 1

            # Test that result has correct type
            @test result_dkmer isa typeof(dkmer)

            # Test that original is unchanged (immutability)
            @test length(dkmer) == length(result_dkmer) + 1
        end
    end

    # Test specific sequences to verify correctness
    @test pop(dmer"TAGCA"d) == dmer"TAGC"d
    @test pop(dmer"AUGCU"r) == dmer"AUGC"r
    @test pop_first(dmer"TAGCA"d) == dmer"AGCA"d
    @test pop_first(dmer"AUGCU"r) == dmer"UGCU"r

    # Test chaining operations
    @test pop(pop(dmer"TAGCA"d)) == dmer"TAG"d
    @test pop_first(pop_first(dmer"TAGCA"d)) == dmer"GCA"d

    # Test popping down to empty
    @test pop(dmer"T"d) == dmer""d
    @test pop_first(dmer"A"r) == dmer""r

    # Test error when popping empty kmer
    @test_throws BoundsError pop(dmer""d)
    @test_throws BoundsError pop_first(dmer""d)
    @test_throws BoundsError pop(DynamicAAKmer{UInt64}(aa""))
    @test_throws BoundsError pop_first(DynamicAAKmer{UInt128}(aa""))

    # Test with larger bit integers
    d512 = DynamicDNAKmer{UInt512}(dna"TAGCA")
    @test pop(d512) == DynamicDNAKmer{UInt512}(dna"TAGC")
    @test pop_first(d512) == DynamicDNAKmer{UInt512}(dna"AGCA")
end

@testset "setindex" begin
    # Basic functionality
    d = dmer"TAGCTGA"d
    @test Base.setindex(d, DNA_C, 1) == dmer"CAGCTGA"d
    @test Base.setindex(d, DNA_A, 4) == dmer"TAGATGA"d
    @test Base.setindex(d, DNA_T, 7) == dmer"TAGCTGT"d
    @test d == dmer"TAGCTGA"d  # Original unchanged

    # Different alphabets and backing types
    @test Base.setindex(dmer"AUGC"r, RNA_G, 2) == dmer"AGGC"r
    @test Base.setindex(dmer"KWOP"a, AA_L, 3) == dmer"KWLP"a
    @test Base.setindex(DynamicDNAKmer{UInt32}(dna"ATGC"), DNA_G, 2) == DynamicDNAKmer{UInt32}(dna"AGGC")

    # Type conversion
    @test Base.setindex(dmer"TAG"d, 'C', 2) == dmer"TCG"d

    # Bounds checking
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 0)
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 5)
    @test_throws BoundsError Base.setindex(dmer""d, DNA_A, 1)

    # Test with larger bit integers
    d256 = DynamicDNAKmer{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_G, 3) == DynamicDNAKmer{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_C, 1) == DynamicDNAKmer{UInt256}(dna"CAGCAT")
end

@testset "capacity" begin
    # Create a zero-BPS alphabet for testing
    struct ZeroBPSAlphabet <: Alphabet end
    Base.eltype(::Type{ZeroBPSAlphabet}) = DNA
    BioSequences.BitsPerSymbol(::ZeroBPSAlphabet) = BioSequences.BitsPerSymbol{0}()

    # Test zero BPS: capacity should be clamped to typemax(Int)
    @test capacity(DynamicKmer{ZeroBPSAlphabet, UInt8}) == clamp(typemax(UInt8), Int)
    @test capacity(DynamicKmer{ZeroBPSAlphabet, UInt32}) == clamp(typemax(UInt32), Int)
    @test capacity(DynamicKmer{ZeroBPSAlphabet, UInt128}) == typemax(Int)

    # Test non-zero BPS: capacity should be in range 0:div(8 * sizeof(U), B)
    for (A, bps) in [(DNAAlphabet{2}, 2), (DNAAlphabet{4}, 4), (AminoAcidAlphabet, 8)]
        for U in [UInt8, UInt32, UInt64]
            cap = capacity(DynamicKmer{A, U})
            max_possible = div(8 * sizeof(U), bps)
            @test 0 <= cap <= max_possible
        end
    end

    # Test that larger backing type gives larger or equal capacity
    @test capacity(DynamicDNAKmer{UInt32}) <= capacity(DynamicDNAKmer{UInt64})
    @test capacity(DynamicAAKmer{UInt32}) <= capacity(DynamicAAKmer{UInt64})
    @test capacity(DynamicDNAKmer{UInt64}) <= capacity(DynamicDNAKmer{UInt256})
    @test capacity(DynamicAAKmer{UInt128}) <= capacity(DynamicAAKmer{UInt512})
end

@testset "translate" begin
    # Basic translation - 2-bit alphabets
    @testset "Basic 2-bit" begin
        # DNA 2-bit with different backing types
        @test translate(dmer"ATGTAA"d) == dmer"M*"a
        @test translate(DynamicDNAKmer{UInt32}(dna"ATGTAA")) == dmer"M*"a
        @test translate(DynamicDNAKmer{UInt64}(dna"ATGTAA")) == dmer"M*"a

        # RNA 2-bit with different backing types
        @test translate(dmer"AUGUAA"r) == dmer"M*"a
        @test translate(DynamicRNAKmer{UInt32}(rna"AUGUAA")) == dmer"M*"a
        @test translate(DynamicRNAKmer{UInt64}(rna"AUGUAA")) == dmer"M*"a

        # Longer sequences
        @test translate(dmer"TCTACACCCTAG"d) == dmer"STP*"a
        @test translate(dmer"UCUACACCCUAG"r) == dmer"STP*"a
    end

    # Basic translation - 4-bit alphabets
    @testset "Basic 4-bit" begin
        # DNA 4-bit with different backing types
        d = DynamicKmer{DNAAlphabet{4}, UInt64}(dna"TGGCCCGATTGA")
        @test translate(d) == dmer"WPD*"a

        # UInt128 works with 4-bit since capacity is smaller
        d128 = DynamicKmer{DNAAlphabet{4}, UInt128}(dna"ATGTAA")
        @test translate(d128) == dmer"M*"a

        # RNA 4-bit with different backing types
        r = DynamicKmer{RNAAlphabet{4}, UInt64}(rna"UGGCCCGAUUGA")
        @test translate(r) == dmer"WPD*"a

        r32 = DynamicKmer{RNAAlphabet{4}, UInt32}(rna"AUGUAA")
        @test translate(r32) == dmer"M*"a

        r128 = DynamicKmer{RNAAlphabet{4}, UInt128}(rna"AUGUAA")
        @test translate(r128) == dmer"M*"a
    end

    # Different genetic codes
    @testset "Different genetic codes" begin
        # Vertebrate mitochondrial code: AGA and AGG are stop codons
        vert_mito = BioSequences.vertebrate_mitochondrial_genetic_code
        @test translate(dmer"ATGAGA"d; code = vert_mito) == dmer"M*"a  # AGA is stop

        # Standard code: AGA and AGG code for R (Arginine)
        @test translate(dmer"ATGAGA"d) == dmer"MR"a
        @test translate(dmer"ATGAGG"d) == dmer"MR"a

        # Test with different backing types
        @test translate(DynamicDNAKmer{UInt64}(dna"ATGAGA"); code = vert_mito) == dmer"M*"a
        @test translate(DynamicRNAKmer{UInt32}(rna"AUGAGA"); code = vert_mito) == dmer"M*"a

        # Test with 4-bit alphabets
        @test translate(DynamicKmer{DNAAlphabet{4}, UInt64}(dna"ATGAGA"); code = vert_mito) == dmer"M*"a
    end

    # alternative_start flag
    @testset "alternative_start" begin
        # Without alternative_start: TTG codes for L (Leucine)
        @test translate(dmer"TTGCCC"d; alternative_start = false) == dmer"LP"a
        @test translate(dmer"UUGCCC"r; alternative_start = false) == dmer"LP"a

        # With alternative_start: first codon becomes M regardless
        @test translate(dmer"TTGCCC"d; alternative_start = true) == dmer"MP"a
        @test translate(dmer"UUGCCC"r; alternative_start = true) == dmer"MP"a
    end

    # Ambiguous codons (only for 4-bit alphabets)
    @testset "Ambiguous codons" begin
        # With allow_ambiguous_codons=true (default), ambiguous codons translate
        d_ambig = DynamicKmer{DNAAlphabet{4}, UInt128}("AAAACWGCSWTARACADA")
        @test translate(d_ambig) == dmer"KTAJBX"a

        # With allow_ambiguous_codons=false, ambiguous codons throw
        @test_throws Exception translate(d_ambig; allow_ambiguous_codons = false)

        # Test various ambiguous nucleotides
        # W = A or T, so TWG could be AAG (K) or TAG (*)
        # With allow_ambiguous, should give X (ambiguous)
        d_w = DynamicKmer{DNAAlphabet{4}, UInt64}(dna"ATGTWG")
        result_w = translate(d_w; allow_ambiguous_codons = true)
        @test length(result_w) == 2  # M and something
    end

    # Error: Length not divisible by 3
    @testset "Length not divisible by 3" begin
        @test_throws ArgumentError translate(dmer"A"d)
        @test_throws ArgumentError translate(dmer"UG"r)
        @test_throws ArgumentError translate(dmer"CUGUAGUUGUCGC"r)
        @test_throws ArgumentError translate(DynamicKmer{RNAAlphabet{4}, UInt32}(rna"AUGC"))
    end

    # Error: Sequences with gaps (only for 4-bit alphabets)
    @testset "Sequences with gaps" begin
        @test_throws Exception translate(DynamicKmer{RNAAlphabet{4}, UInt64}(rna"-UGAUG"))
        @test_throws Exception translate(DynamicKmer{DNAAlphabet{4}, UInt64}(dna"AT-ATG"))
        @test_throws Exception translate(DynamicKmer{RNAAlphabet{4}, UInt64}(rna"AUGAU-"))
        @test_throws Exception translate(DynamicKmer{DNAAlphabet{4}, UInt64}(dna"A--"))
    end

    # Error: Input type capacity too large for output
    @testset "Capacity overflow" begin
        # UInt128 with 2-bit alphabet has capacity ~63, which translates to 21 AA
        # needing 168 bits, exceeding UInt128's 128 bits
        # This should error at translation time
        @test_throws ErrorException translate(DynamicDNAKmer{UInt128}(dna"ATG"))
        @test_throws ErrorException translate(DynamicRNAKmer{UInt128}(rna"AUG"))

        # Even empty sequences should error due to type-based capacity check
        @test_throws ErrorException translate(DynamicDNAKmer{UInt128}(dna""))
        @test_throws ErrorException translate(DynamicRNAKmer{UInt128}(rna""))
    end

    # Edge cases
    @testset "Edge cases" begin
        # Empty sequence (length 0 is divisible by 3, produces 0 AA)
        @test translate(dmer""d) == dmer""a
        @test translate(dmer""r) == dmer""a
        @test translate(DynamicDNAKmer{UInt32}(dna"")) == dmer""a
        @test translate(DynamicKmer{DNAAlphabet{4}, UInt64}(dna"")) == dmer""a

        # Very long sequence (near capacity) - 2-bit
        # For UInt64 with 2-bit DNA, capacity is 31, so 30 nt = 10 AA
        long_dna = dmer"ATGATGATGATGATGATGATGATGATG"d
        @test length(translate(long_dna)) == 9
        @test LongSequence(translate(long_dna)) == aa"MMMMMMMMM"
    end
end

@testset "Misc" begin
    d = DynamicAAKmer{UInt32}("WPK")
    @test only([d]') === d
end
