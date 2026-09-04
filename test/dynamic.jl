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
                m = Oligo{typeof(Alphabet(s)), UInt64}(i)
                @test length(m) == length(i)
                @test string(m) == string(s)
            end
        end

        # Test macro equivalence with constructor
        @test dmer"ATGTCGTTAGT"d == DNAOligo{UInt64}(dna"ATGTCGTTAGT")
        @test dmer"AUGUCGAGUGUAUGC"r == RNAOligo{UInt64}(rna"AUGUCGAGUGUAUGC")
        @test dmer"KWOP"a == AAOligo{UInt64}(aa"KWOP")
        @test_throws Exception eval(:(dmer"ATCGATAG"k))

        d = DNAOligo{UInt}(dna"ATGTCGTTAGT")
        @test RNAOligo(d) == d
        @test RNAOligo{UInt64}(d) == d

        # From a large kmer
        m = mer"TAGTGCTGTAGTAGTGCTGTATGATGTCTGCATGC"d
        dm = DNAOligo{UInt256}(m)
        @test LongSequence(m) == LongSequence(dm)

        # Using large bitintegers - 240 2-bit nt = 480 bits should fit in
        # 512 bits with plenty room for the runtime length
        seq = randdnaseq(240)
        dm = DNAOligo{UInt512}(seq)
        @test string(seq) == string(dm)
    end

    @testset "Default backing constructors" begin
        for (constructor, source, expected_type) in [
                (Oligo{DNAAlphabet{2}}, "TAG", DNAOligo{UInt64}),
                (DNAOligo, "TAG", DNAOligo{UInt64}),
                (RNAOligo, rna"AUG", RNAOligo{UInt64}),
                (AAOligo, "KWOP", AAOligo{UInt64}),
                (Oligo{DNAAlphabet{4}}, "TAGN", Oligo{DNAAlphabet{4}, UInt64}),
            ]
            oligomer = @inferred constructor(source)
            @test typeof(oligomer) === expected_type
        end

        narrow = DNAOligo{UInt8}("TAG")
        @test @inferred(DNAOligo(narrow)) === DNAOligo{UInt64}(narrow)
        @test @inferred(RNAOligo(narrow)) === RNAOligo{UInt64}(narrow)
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
                    dst = Oligo{A{dstB}, UInt64}(src)

                    @test string(src) == string(dst)
                end
            end
        end
    end

    @testset "Four to two bit alphabet" begin
        for s in [dna"TAGCTGAC", dna"ATGCTA", dna""]
            for A in [DNAAlphabet{2}, RNAAlphabet{2}]
                m = Oligo{A, UInt64}(LongSequence{DNAAlphabet{4}}(s))
                @test string(m) == string(LongSequence{A}(s))
            end
        end
    end

    @testset "Generic alphabet" begin
        for s in ["HE", "", "中Å!"]
            m = Oligo{CharAlphabet, UInt512}(s)
            @test length(m) == length(s)
            @test string(m) == s
        end
    end

    @testset "From Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(s)
            dkmer = Oligo{typeof(Alphabet(s)), UInt64}(kmer)
            @test string(dkmer) == string(s)
            @test length(dkmer) == length(kmer)
        end
    end

    @testset "From Kmer into wider backing types" begin
        for U in [UInt128, UInt256, UInt512]
            # One symbol, a partially filled word, and a full source word.
            for s in [dna"T", dna"TAG", dna"T"^32]
                kmer = DNAKmer{length(s)}(s)
                oligomer = DNAOligo{U}(kmer)
                @test string(oligomer) == string(kmer)
                @test DNAKmer{length(s)}(oligomer) == kmer
            end

            # Keep a multiword source case covered too.
            s = dna"T"^33
            kmer = DNAKmer{33}(s)
            oligomer = DNAOligo{U}(kmer)
            @test string(oligomer) == string(kmer)
            @test DNAKmer{33}(oligomer) == kmer
        end
    end

    @testset "From Oligo" begin
        # Same backing type, same alphabet
        s = dna"TAGCTGA"
        d1 = dmer"TAGCTGA"d
        d2 = DNAOligo{UInt64}(d1)
        @test d2 == d1
        @test d2 === d1  # Should be identical
        @test length(d2) == length(d1)

        # Different backing type, same alphabet - widening
        d32 = DNAOligo{UInt32}(dna"TAGC")
        d64 = DNAOligo{UInt64}(d32)
        @test_throws MethodError d64 == d32
        @test_throws MethodError d32 == d64
        @test length(d64) == length(d32)
        @test string(d64) == "TAGC"

        # Different backing type, same alphabet - narrowing
        d128 = DNAOligo{UInt128}(dna"TAGC")
        d32_narrow = DNAOligo{UInt32}(d128)
        @test_throws MethodError d32_narrow == d128
        @test_throws MethodError d128 == d32_narrow
        @test length(d32_narrow) == length(d128)
        @test string(d32_narrow) == "TAGC"

        # Same alphabet family (DNA/RNA), different alphabets
        d_dna = dmer"ATGTCGTTAGT"d
        d_rna = RNAOligo{UInt64}(d_dna)
        @test d_rna == d_dna
        @test length(d_rna) == length(d_dna)

        # Different backing type AND different alphabet
        d_dna32 = DNAOligo{UInt32}(dna"TAGC")
        d_rna64 = RNAOligo{UInt64}(d_dna32)
        @test_throws MethodError d_rna64 == d_dna32
        @test_throws MethodError d_dna32 == d_rna64
        @test length(d_rna64) == length(d_dna32)

        # Test with amino acids
        aa_d64 = dmer"KWOP"a
        aa_d128 = AAOligo{UInt128}(aa_d64)
        @test_throws MethodError aa_d128 == aa_d64
        @test length(aa_d128) == length(aa_d64)
        aa_d256 = AAOligo{UInt256}(aa_d64)
        @test_throws MethodError aa_d256 == aa_d64
        @test length(aa_d256) == length(aa_d64)

        # Test with empty kmers
        empty_d32 = DNAOligo{UInt32}(dna"")
        empty_d64 = DNAOligo{UInt64}(empty_d32)
        @test_throws MethodError empty_d64 == empty_d32
        @test length(empty_d64) == 0

        # Test with 4-bit alphabets
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCN")
        d_4bit_32 = Oligo{DNAAlphabet{4}, UInt32}(s_4bit)
        d_4bit_64 = Oligo{DNAAlphabet{4}, UInt64}(d_4bit_32)
        @test_throws MethodError d_4bit_64 == d_4bit_32
        @test length(d_4bit_64) == length(d_4bit_32)

        # Narrowing must compare symbol capacities, not coding bits with symbols.
        for (source, target, fitting, overflowing) in [
                (DNAOligo{UInt16}, DNAOligo{UInt8}, dna"TAG", dna"TAGC"),
                (RNAOligo{UInt16}, RNAOligo{UInt8}, rna"UAG", rna"UAGC"),
            ]
            original = source(fitting)
            narrowed = target(original)
            @test string(narrowed) == string(original)
            @test length(narrowed) == length(original)
            @test_throws MethodError narrowed == original
            @test_throws Exception target(source(overflowing))
        end

        # AAOligo{UInt8} has capacity zero, exercising the smallest target
        # backing type too.
        empty_aa = AAOligo{UInt16}(aa"")
        @test_throws MethodError AAOligo{UInt8}(empty_aa) == empty_aa
        @test_throws Exception AAOligo{UInt8}(AAOligo{UInt16}(aa"K"))
    end

    @testset "Widening" begin
        for (source, target) in [
                (DNAOligo{UInt8}, DNAOligo{UInt16}),
                (RNAOligo{UInt16}, RNAOligo{UInt32}),
                (AAOligo{UInt32}, AAOligo{UInt64}),
                (Oligo{DNAAlphabet{4}, UInt64}, Oligo{DNAAlphabet{4}, UInt128}),
                (DNAOligo{UInt128}, DNAOligo{UInt256}),
                (RNAOligo{UInt256}, RNAOligo{UInt512}),
                (AAOligo{UInt512}, AAOligo{UInt1024}),
            ]
            @test widen(source) === target
        end

        for mer in [
                DNAOligo{UInt8}(dna"TAG"),
                RNAOligo{UInt32}(rna"UAGC"),
                AAOligo{UInt128}(aa"KWOP"),
                Oligo{DNAAlphabet{4}, UInt256}(dna"TAGCN"),
            ]
            widened = widen(mer)
            @test typeof(widened) === widen(typeof(mer))
            @test widened === widen(typeof(mer))(mer)
            @test string(widened) == string(mer)
            @test length(widened) == length(mer)
        end

        @test_throws MethodError widen(DNAOligo{UInt1024})
    end

    @testset "To Kmer" begin
        for s in [dna"TAGCTA", rna"UGCUGA", aa"PLKWM"]
            dkmer = Oligo{typeof(Alphabet(s)), UInt64}(s)
            kmer = Kmer{typeof(Alphabet(s)), length(s)}(dkmer)
            @test length(kmer) == length(dkmer)
            @test string(dkmer) == string(kmer)
            @test_throws MethodError kmer == dkmer
            @test_throws MethodError dkmer == kmer
        end
    end

    @testset "To Kmer with N=2 (>64 coding bits)" begin
        # For 2-bit DNA: need >32 bases for >64 coding bits
        s_dna = dna"TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"  # 33 bases = 66 bits
        dkmer_dna = Oligo{DNAAlphabet{2}, UInt256}(s_dna)
        kmer_dna = Kmer{DNAAlphabet{2}, 33}(dkmer_dna)
        @test length(kmer_dna) == length(dkmer_dna)
        @test string(dkmer_dna) == string(kmer_dna)

        # For 2-bit RNA
        s_rna = rna"UGCUGAUGCUGAUGCUGAUGCUGAUGCUGAUGA"  # 33 bases = 66 bits
        dkmer_rna = Oligo{RNAAlphabet{2}, UInt256}(s_rna)
        kmer_rna = Kmer{RNAAlphabet{2}, 33}(dkmer_rna)
        @test length(kmer_rna) == length(dkmer_rna)
        @test string(dkmer_rna) == string(kmer_rna)

        # For 8-bit amino acids: need >8 bases for >64 coding bits
        s_aa = aa"KWOPPLKWM"  # 9 bases = 72 bits
        dkmer_aa = Oligo{AminoAcidAlphabet, UInt256}(s_aa)
        kmer_aa = Kmer{AminoAcidAlphabet, 9}(dkmer_aa)
        @test length(kmer_aa) == length(dkmer_aa)
        @test string(dkmer_aa) == string(kmer_aa)

        # Test error on length mismatch
        dkmer = dmer"TAG"d
        @test_throws Exception Kmer{DNAAlphabet{2}, 5}(dkmer)
    end

    @testset "Capacity limits" begin
        # Test that exceeding capacity throws
        @test_throws Exception DNAOligo{UInt32}(dna"T"^30)
        @test_throws Exception AAOligo{UInt32}(aa"A"^8)
    end
end

@testset "Unsafe runtime-length extraction" begin
    unsafe_extract = Kmers.unsafe_extract

    copy_source = LongDNA{2}("CCTAGCAAT")
    @test @inferred(
        unsafe_extract(Kmers.Copyable(), DNAOligo{UInt32}, copy_source, 3, 4)
    ) == DNAOligo{UInt32}("TAGC")
    @test @inferred(
        unsafe_extract(Kmers.Copyable(), DNAOligo{UInt32}, copy_source, 10, 0)
    ) === empty(DNAOligo{UInt32})

    two_bit_source = LongRNA{2}("CCUAGCAAU")
    four_bit_type = Oligo{DNAAlphabet{4}, UInt32}
    @test @inferred(
        unsafe_extract(Kmers.TwoToFour(), four_bit_type, two_bit_source, 3, 4)
    ) == four_bit_type("TAGC")
    @test @inferred(
        unsafe_extract(Kmers.TwoToFour(), four_bit_type, two_bit_source, 10, 0)
    ) === empty(four_bit_type)

    four_bit_source = LongDNA{4}("CCTAGCAAN")
    @test @inferred(
        unsafe_extract(Kmers.FourToTwo(), DNAOligo{UInt32}, four_bit_source, 3, 4)
    ) == DNAOligo{UInt32}("TAGC")
    @test @inferred(
        unsafe_extract(Kmers.FourToTwo(), DNAOligo{UInt32}, four_bit_source, 10, 0)
    ) === empty(DNAOligo{UInt32})
    @test_throws BioSequences.EncodeError unsafe_extract(
        Kmers.FourToTwo(), DNAOligo{UInt32}, four_bit_source, 7, 3
    )

    ascii_source = codeunits("XXTAGCXP")
    @test @inferred(
        unsafe_extract(Kmers.AsciiEncode(), DNAOligo{UInt32}, ascii_source, 3, 4)
    ) == DNAOligo{UInt32}("TAGC")
    @test @inferred(
        unsafe_extract(Kmers.AsciiEncode(), DNAOligo{UInt32}, ascii_source, 10, 0)
    ) === empty(DNAOligo{UInt32})
    @test_throws BioSequences.EncodeError unsafe_extract(
        Kmers.AsciiEncode(), DNAOligo{UInt32}, ascii_source, 8, 1
    )

    generic_source = LongSequence{CharAlphabet}("中Å!人大网")
    generic_type = Oligo{CharAlphabet, UInt256}
    @test @inferred(
        unsafe_extract(Kmers.GenericRecoding(), generic_type, generic_source, 2, 3)
    ) == generic_type("Å!人")
    @test @inferred(
        unsafe_extract(Kmers.GenericRecoding(), generic_type, generic_source, 7, 0)
    ) === empty(generic_type)
end

@testset "Packed unsafe runtime-length extraction" begin
    unsafe_extract = Kmers.unsafe_extract

    function compare_packed_extract(recoding, T, source, from, len)
        @test unsafe_extract(recoding, T, source, from, len) ==
            Kmers.extract_elements(recoding, T, source, from, len)
    end

    source_2bit = LongDNA{2}(repeat("ACGT", 100))
    source_4bit = LongDNA{4}(repeat("ACGT", 100))
    for U in (UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512),
            (recoding, A, source) in (
                (Kmers.Copyable(), DNAAlphabet{2}, source_2bit),
                (Kmers.TwoToFour(), DNAAlphabet{4}, source_2bit),
                (Kmers.FourToTwo(), DNAAlphabet{2}, source_4bit),
            )
        T = Oligo{A, U}
        for len in (0, 1, min(3, capacity(T)), capacity(T)), from in (1, 2, 3, 31, 32)
            (len <= capacity(T) && from + len - 1 <= length(source)) || continue
            compare_packed_extract(recoding, T, source, from, len)
        end
    end

    # Exercise copyable packed extraction with non-nucleic and one-bit alphabets.
    for (A, source) in (
                (AminoAcidAlphabet, LongAA(repeat("ARNDCQEGHILKMFPSTWYV", 20))),
                (
                    OneBPSAlphabet,
                    LongSequence{OneBPSAlphabet}(repeat([DNA_A, DNA_C], 200)),
                ),
            ), U in (UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512)
        T = Oligo{A, U}
        for len in (0, 1, min(3, capacity(T)), capacity(T)), from in (1, 2, 31, 32)
            from + len - 1 <= length(source) || continue
            compare_packed_extract(Kmers.Copyable(), T, source, from, len)
        end
    end

    # DNA and RNA with the same bit width can take the direct-copy path.
    compare_packed_extract(
        Kmers.Copyable(),
        RNAOligo{UInt128},
        source_2bit,
        2,
        33,
    )

    # LongSubSeq preserves its data offset, including across physical-word boundaries.
    source_2bit_view = view(source_2bit, 2:300)
    source_4bit_view = view(source_4bit, 2:300)
    for len in (0, 1, 16, 17, 32, 33), from in (1, 2, 31, 32)
        if len <= capacity(DNAOligo{UInt128})
            compare_packed_extract(
                Kmers.Copyable(),
                DNAOligo{UInt128},
                source_2bit_view,
                from,
                len,
            )
            compare_packed_extract(
                Kmers.FourToTwo(),
                DNAOligo{UInt128},
                source_4bit_view,
                from,
                len,
            )
        end
        if len <= capacity(Oligo{DNAAlphabet{4}, UInt128})
            compare_packed_extract(
                Kmers.TwoToFour(),
                Oligo{DNAAlphabet{4}, UInt128},
                source_2bit_view,
                from,
                len,
            )
        end
    end

    function extract_error(f)
        try
            f()
        catch error
            return error
        end
        error("Expected unsafe_extract to throw")
    end

    # The packed certainty check is range-local and reports the first uncertain
    # or gap symbol in the rejected chunk.
    for (position, symbol) in ((1, 'N'), (17, '-'), (33, 'N'))
        chars = fill('A', 33)
        chars[position] = symbol
        source = LongDNA{4}(String(chars))
        error = extract_error() do
            unsafe_extract(Kmers.FourToTwo(), DNAOligo{UInt128}, source, 1, 33)
        end
        @test error isa BioSequences.EncodeError
        @test error.val == source[position]
    end
    outside_range = LongDNA{4}("N" * repeat("ACGT", 8) * "N")
    @test unsafe_extract(Kmers.FourToTwo(), DNAOligo{UInt128}, outside_range, 2, 32) ==
        Kmers.extract_elements(Kmers.FourToTwo(), DNAOligo{UInt128}, outside_range, 2, 32)

    function packed_extract_allocations(source_2bit, source_4bit)
        two_bit_type = DNAOligo{UInt128}
        four_bit_type = Oligo{DNAAlphabet{4}, UInt128}
        unsafe_extract(Kmers.Copyable(), two_bit_type, source_2bit, 2, 33)
        unsafe_extract(Kmers.TwoToFour(), four_bit_type, source_2bit, 2, 17)
        unsafe_extract(Kmers.FourToTwo(), two_bit_type, source_4bit, 2, 33)
        return (
            @allocated(unsafe_extract(Kmers.Copyable(), two_bit_type, source_2bit, 2, 33)),
            @allocated(unsafe_extract(Kmers.TwoToFour(), four_bit_type, source_2bit, 2, 17)),
            @allocated(unsafe_extract(Kmers.FourToTwo(), two_bit_type, source_4bit, 2, 33)),
        )
    end
    @test all(iszero, packed_extract_allocations(source_2bit, source_4bit))
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
        @test m[9:8] == dmer""d

        @test_throws BoundsError m[0:3]
        @test_throws BoundsError m[6:9]
        @test_throws BoundsError m[0:-1]
        @test_throws BoundsError m[10:9]

        empty_mer = dmer""d
        @test @inferred(empty_mer[1:0]) === empty_mer
        @test_throws BoundsError empty_mer[0:-1]
        @test_throws BoundsError empty_mer[2:1]
    end

    @testset "Integer-vector indexing" begin
        m = dmer"TAGCTGAC"d
        @test @inferred(m[[1, 3, 5]]) == dmer"TGT"d
        @test @inferred(m[[8, 2, 8, 1]]) == dmer"CACT"d
        @test @inferred(m[Int[]]) === empty(typeof(m))

        repeated = DNAOligo{UInt8}(dna"TAG")
        @test @inferred(repeated[[3, 3, 1]]) == DNAOligo{UInt8}(dna"GGT")
        @test_throws BoundsError repeated[[1, 1, 1, 1]]
        @test_throws BoundsError (@inbounds repeated[[1, 1, 1, 1]])

        for (x, indices, expected) in Any[
                (Oligo{RNAAlphabet{4}, UInt64}(rna"AUGN"), [4, 2, 4], rna"NUN"),
                (AAOligo{UInt128}(aa"KWOP"), [4, 1, 1], aa"PKK"),
                (Oligo{CharAlphabet, UInt512}("中Å!"), [3, 1, 2], "!中Å"),
            ]
            @test string(x[indices]) == string(expected)
            @test typeof(x[indices]) === typeof(x)
        end

        @test_throws BoundsError m[[2, 9]]
        @test_throws BoundsError m[[0, 1]]
        @test_throws BoundsError m[[-1]]
    end

    @testset "Logical indexing" begin
        m = dmer"TAGCTGAC"d
        @test @inferred(m[Bool[1, 0, 1, 0, 0, 1, 0, 1]]) == dmer"TGGC"d
        @test @inferred(m[trues(length(m))]) === m
        @test @inferred(m[falses(length(m))]) === empty(typeof(m))

        empty_mer = AAOligo{UInt64}(aa"")
        @test @inferred(empty_mer[Bool[]]) === empty_mer

        @test_throws BoundsError m[trues(length(m) + 1)]
        @test_throws BoundsError m[falses(length(m) - 1)]
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

    @testset "Comparison rejects different backing types" begin
        s1 = dna"TAG"
        s2 = dna"TAC"
        s3 = dna"TAGA"

        m32_1 = DNAOligo{UInt32}(s1)
        m64_1 = DNAOligo{UInt64}(s1)
        m32_2 = DNAOligo{UInt32}(s2)
        m64_2 = DNAOligo{UInt64}(s2)

        for (a, b) in [(m32_1, m64_1), (m32_2, m64_2)]
            @test_throws MethodError a == b
            @test_throws MethodError b == a
            @test_throws MethodError isequal(a, b)
            @test_throws MethodError a < b
            @test_throws MethodError cmp(a, b)
        end

        # Different lengths do not make the backing types comparable.
        m32_3 = DNAOligo{UInt32}(s3)
        @test m32_1 < m32_3
        @test cmp(m32_1, m32_3) < 0
        @test_throws MethodError m64_1 < m32_3
        @test_throws MethodError cmp(m64_1, m32_3)

        # Also for AA types
        s1 = aa"KPRCRLF"
        s2 = aa"KPRCRLFAAA"

        m64_1 = AAOligo{UInt64}(s1)
        m128_1 = AAOligo{UInt128}(s1)
        m128_2 = AAOligo{UInt128}(s2)
        m256_1 = AAOligo{UInt256}(s1)

        @test_throws MethodError m64_1 == m128_1
        @test_throws MethodError cmp(m64_1, m128_1)
        @test m128_1 < m128_2
        @test_throws MethodError m64_1 == m256_1
        @test_throws MethodError cmp(m256_1, m128_1)
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
            m = DNAOligo{UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m == m2
        end

        # Error on exceeding capacity
        @test_throws Exception from_integer(DNAOligo{UInt32}, UInt32(0), 30)

        # A zero-length value must discard every input bit, including inputs
        # which would otherwise wrap a full-width shift to a zero-bit shift.
        for U in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512]
            T = DNAOligo{U}
            for value in [zero(U), one(U), typemax(U)]
                oligomer = from_integer(T, value, 0)
                @test length(oligomer) == 0
                @test isempty(oligomer)
                @test oligomer == empty(T)
                @test oligomer.x == zero(U)
            end
        end
    end

    @testset "Round-trip conversion" begin
        for s in [dna"TAGCTGA", rna"UGCUGA", aa"PLKWM"]
            m = Oligo{typeof(Alphabet(s)), UInt64}(s)
            u = as_integer(m)
            m2 = from_integer(typeof(m), u, length(m))
            @test m === m2
        end

        # Use a short sequence whose coding bits fit in every input width, and
        # vary the input and target widths independently.
        for Target in [UInt8, UInt16, UInt32, UInt64, UInt128],
                Input in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256]
            m = DNAOligo{Target}(dna"TAG")
            input = as_integer(m) % Input
            @test from_integer(typeof(m), input, length(m)) === m
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
        @test Kmers.fx_hash(m1, UInt(0)) == Kmers.fx_hash(m2, UInt(0))

        # Different kmers should produce different fx_hash values
        @test Kmers.fx_hash(m1, UInt(0)) != Kmers.fx_hash(m3, UInt(0))

        # Different seeds should produce different hashes
        @test Kmers.fx_hash(m1, UInt(123)) != Kmers.fx_hash(m1, UInt(456))

        # Length must be part of hash: TCA and TC have identical coding bits
        # (since A encodes to 00, which looks like padding), but different lengths.
        # They must hash differently despite having the same as_integer representation.
        m1 = dmer"TCA"d
        m2 = dmer"TC"d
        @test Kmers.fx_hash(m1, UInt(0)) != Kmers.fx_hash(m2, UInt(0))
    end
end

@testset "Biological operations" begin
    @testset "Reverse" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test reverse(m) == DNAOligo{UInt64}(reverse(s))

        # Test empty sequence
        @test reverse(dmer""d) == dmer""d

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse(m_rna) == RNAOligo{UInt64}(reverse(s_rna))

        # Test with 4-bit DNA
        s_4bit = dna"TAGCTGA"
        m_4bit = Oligo{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse(m_4bit) == Oligo{DNAAlphabet{4}, UInt64}(reverse(s_4bit))

        # Test with amino acids
        s_aa = aa"KWOP"
        m_aa = dmer"KWOP"a
        @test reverse(m_aa) == AAOligo{UInt64}(reverse(s_aa))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGC"
        m_large = DNAOligo{UInt256}(s_large)
        @test reverse(m_large) == DNAOligo{UInt256}(reverse(s_large))
    end

    @testset "Complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test complement(m) == DNAOligo{UInt64}(complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test complement(m_rna) == RNAOligo{UInt64}(complement(s_rna))

        # Test with 4-bit DNA (includes ambiguous bases)
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGAW")  # W = A or T
        m_4bit = Oligo{DNAAlphabet{4}, UInt64}(s_4bit)
        @test complement(m_4bit) == Oligo{DNAAlphabet{4}, UInt64}(complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGAW")  # W = A or U
        m_rna_4bit = Oligo{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test complement(m_rna_4bit) == Oligo{RNAAlphabet{4}, UInt64}(complement(s_rna_4bit))

        # Test generic fallback with a nonstandard nucleic-acid alphabet.
        generic = Oligo{GenericNucAlphabet, UInt64}(dna"TAGC")
        @test complement(generic) == Oligo{GenericNucAlphabet, UInt64}(dna"ATCG")
        @test complement(empty(generic)) == empty(generic)
    end

    @testset "Reverse complement" begin
        # Test with 2-bit DNA
        s = dna"TAGCTGA"
        m = dmer"TAGCTGA"d
        @test reverse_complement(m) == DNAOligo{UInt64}(reverse_complement(s))

        # Test with 2-bit RNA
        s_rna = rna"UAGCUGA"
        m_rna = dmer"UAGCUGA"r
        @test reverse_complement(m_rna) == RNAOligo{UInt64}(reverse_complement(s_rna))

        # Test with 4-bit DNA
        s_4bit = LongSequence{DNAAlphabet{4}}(dna"TAGCTGA")
        m_4bit = Oligo{DNAAlphabet{4}, UInt64}(s_4bit)
        @test reverse_complement(m_4bit) == Oligo{DNAAlphabet{4}, UInt64}(reverse_complement(s_4bit))

        # Test with 4-bit RNA
        s_rna_4bit = LongSequence{RNAAlphabet{4}}(rna"UAGCUGA")
        m_rna_4bit = Oligo{RNAAlphabet{4}, UInt64}(s_rna_4bit)
        @test reverse_complement(m_rna_4bit) == Oligo{RNAAlphabet{4}, UInt64}(reverse_complement(s_rna_4bit))

        # Test with larger bit integers
        s_large = dna"TAGCTAGCTAGCTAGCTAGC"
        m_large = DNAOligo{UInt512}(s_large)
        @test reverse_complement(m_large) == DNAOligo{UInt512}(reverse_complement(s_large))

        generic = Oligo{GenericNucAlphabet, UInt64}(dna"TAGC")
        @test reverse_complement(generic) == Oligo{GenericNucAlphabet, UInt64}(dna"GCTA")
        @test reverse_complement(empty(generic)) == empty(generic)
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
        m_large = DNAOligo{UInt256}(dna"TAGCTAGCTAGC")
        m_rc = reverse_complement(m_large)
        @test canonical(m_large) == min(m_large, m_rc)

        generic = Oligo{GenericNucAlphabet, UInt64}(dna"TAGC")
        @test canonical(generic) == reverse_complement(generic)
        @test canonical(empty(generic)) == empty(generic)
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

        m = AAOligo{UInt256}(aa"WLAKWVMARQKW")
        @test count(==(AA_W), m) == 3
        @test count(==(AA_Q), m) == 1
        @test count(==(AA_A), m) == 2
        @test count(==(AA_C), m) == 0

        # Test with even larger integers
        m_large = DNAOligo{UInt512}(dna"TAGCTAGCTAGCTAGC")
        @test count(==(DNA_T), m_large) == 4
        @test count(==(DNA_A), m_large) == 4

        # Test all standard alphabets and backing widths against explicit
        # iteration. This covers zero encodings, absent symbols, empty values,
        # and values that fill their backing type to capacity.
        for A in (
                    DNAAlphabet{2},
                    RNAAlphabet{2},
                    DNAAlphabet{4},
                    RNAAlphabet{4},
                    AminoAcidAlphabet,
                ), U in (UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512)
            T = Oligo{A, U}
            alphabet_symbols = collect(symbols(A()))
            sequence = [
                alphabet_symbols[mod1(i, length(alphabet_symbols))] for
                    i in 1:capacity(T)
            ]
            oligomer = T(sequence)
            @test length(oligomer) == capacity(T)

            for symbol in alphabet_symbols
                @test count(==(symbol), oligomer) == count(==(symbol), collect(oligomer))
                @test count(==(symbol), empty(T)) == 0
            end
        end
    end
end

@testset "shift_encoding" begin
    m = DNAOligo{UInt32}(dna"TAGA")
    enc = UInt32(BioSequences.encode(DNAAlphabet{2}(), DNA_C))
    m2 = Kmers.shift_encoding(m, enc)
    @test m2 == DNAOligo{UInt32}(dna"AGAC")
    @test length(m2) == length(m)
end

@testset "shift and shift_first" begin
    cases = Any[
        (DNAOligo{UInt8}(dna"TAG"), RNA_C),
        (RNAOligo{UInt16}(rna"AUGC"), DNA_T),
        (Oligo{DNAAlphabet{4}, UInt32}(dna"TAGN"), RNA_Y),
        (Oligo{RNAAlphabet{4}, UInt64}(rna"AUGN"), DNA_Y),
        (AAOligo{UInt128}(aa"KWOP"), 'L'),
        (Oligo{CharAlphabet, UInt512}("中Å!"), 'Æ'),
    ]

    for (x, symbol) in cases
        converted = convert(eltype(x), symbol)
        expected_last = [collect(x)[2:end]; converted]
        expected_first = [converted; collect(x)[1:(end - 1)]]

        shifted = @inferred shift(x, symbol)
        shifted_first = @inferred shift_first(x, symbol)
        @test typeof(shifted) === typeof(x)
        @test typeof(shifted_first) === typeof(x)
        @test length(shifted) == length(x)
        @test length(shifted_first) == length(x)
        @test collect(shifted) == expected_last
        @test collect(shifted_first) == expected_first

        U = Kmers.utype(typeof(x))
        encoding = U(BioSequences.encode(Alphabet(x), converted))
        @test @inferred(Kmers.shift_encoding(x, encoding)) === shifted
        @test @inferred(Kmers.shift_first_encoding(x, encoding)) === shifted_first
    end

    for T in (
            DNAOligo{UInt8},
            RNAOligo{UInt64},
            Oligo{DNAAlphabet{4}, UInt128},
            Oligo{RNAAlphabet{4}, UInt256},
            AAOligo{UInt512},
            Oligo{CharAlphabet, UInt512},
        )
        x = empty(T)
        symbol = first(symbols(Alphabet(T)))
        encoding = Kmers.utype(T)(BioSequences.encode(Alphabet(T), symbol))
        @test @inferred(shift(x, symbol)) === x
        @test @inferred(shift_first(x, symbol)) === x
        @test @inferred(Kmers.shift_encoding(x, encoding)) === x
        @test @inferred(Kmers.shift_first_encoding(x, encoding)) === x
    end

    long_sequence = dna"TAGC"^100
    large = DNAOligo{UInt1024}(long_sequence)
    @test @inferred(shift(large, DNA_A)) == DNAOligo{UInt1024}(
        [collect(long_sequence)[2:end]; DNA_A]
    )
    @test @inferred(shift_first(large, DNA_A)) == DNAOligo{UInt1024}(
        [DNA_A; collect(long_sequence)[1:(end - 1)]]
    )
end

@testset "Mixed integer types" begin
    for U in [UInt8, UInt16, UInt32, UInt64, UInt128, UInt256, UInt512]
        s = dna"TAG"
        m = Oligo{DNAAlphabet{2}, U}(s)
        @test string(m) == string(s)
        @test_throws MethodError m == s
        @test_throws MethodError s == m
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
                (DNAOligo{UInt32}(dna"TAG"), DNA_G),
                (RNAOligo{UInt32}(rna"AUG"), RNA_C),
                (AAOligo{UInt64}(aa"KWOP"), AA_L),
                (AAOligo{UInt256}(aa"MWP"), AA_K),
            ]
            # Apply operation to Oligo
            result_dkmer = dkmer_fn(dkmer, symbol)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq, symbol)

            # Test content against the LongSequence reference.
            @test string(result_dkmer) == string(longseq)
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
    d_full = DNAOligo{UInt32}(dna"T"^14)
    @test_throws BoundsError push(d_full, DNA_A)
    @test_throws BoundsError push_first(d_full, DNA_A)

    # For UInt64 with 2-bit DNA: capacity = 29
    d64_full = DNAOligo{UInt64}(dna"A"^29)
    @test_throws BoundsError push(d64_full, DNA_T)
    @test_throws BoundsError push_first(d64_full, DNA_G)

    # For UInt32 with 8-bit AA: capacity = 3
    aa32_full = AAOligo{UInt32}(aa"WPK")
    @test_throws BoundsError push(aa32_full, AA_L)
    @test_throws BoundsError push_first(aa32_full, AA_M)

    # UInt8 has just enough length bits for a full 2-bit oligomer. Incrementing
    # its packed representation first would wrap the length and corrupt it.
    for (T, sequence, symbol) in [
            (DNAOligo{UInt8}, dna"TAG", DNA_A),
            (RNAOligo{UInt8}, rna"UAG", RNA_U),
        ]
        full = T(sequence)
        original = full.x
        @test length(full) == capacity(T)
        @test_throws BoundsError push(full, symbol)
        @test_throws BoundsError push_first(full, symbol)
        @test full.x == original
    end

    # Verify we can push to capacity-1 without error
    d_almost_full = DNAOligo{UInt64}(dna"T"^28)
    @test length(push(d_almost_full, DNA_A)) == 29
    @test length(push_first(d_almost_full, DNA_G)) == 29

    # Test with larger bit integers
    d256 = DNAOligo{UInt256}(dna"TAGC")
    @test push(d256, DNA_G) == DNAOligo{UInt256}(dna"TAGCG")
    @test push_first(d256, DNA_A) == DNAOligo{UInt256}(dna"ATAGC")
end

@testset "pop and pop_first" begin
    for (dkmer_fn, longseq_fn) in Any[[pop, pop!], [pop_first, popfirst!]]
        for dkmer in Any[
                dmer"TAGCA"d,
                dmer"AUGCU"r,
                dmer"T"d,
                dmer"A"r,
                DNAOligo{UInt32}(dna"TAGG"),
                RNAOligo{UInt32}(rna"AUGC"),
                AAOligo{UInt64}(aa"KWOPL"),
                AAOligo{UInt128}(aa"MWPK"),
            ]
            # Apply operation to Oligo
            result_dkmer = dkmer_fn(dkmer)

            # Apply operation to LongSequence for comparison
            longseq = LongSequence{typeof(Alphabet(dkmer))}(dkmer)
            longseq_fn(longseq)

            # Test content against the LongSequence reference.
            @test string(result_dkmer) == string(longseq)
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
    @test_throws BoundsError pop(AAOligo{UInt64}(aa""))
    @test_throws BoundsError pop_first(AAOligo{UInt128}(aa""))

    # Test with larger bit integers
    d512 = DNAOligo{UInt512}(dna"TAGCA")
    @test pop(d512) == DNAOligo{UInt512}(dna"TAGC")
    @test pop_first(d512) == DNAOligo{UInt512}(dna"AGCA")
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
    @test Base.setindex(DNAOligo{UInt32}(dna"ATGC"), DNA_G, 2) == DNAOligo{UInt32}(dna"AGGC")

    # Type conversion
    @test Base.setindex(dmer"TAG"d, 'C', 2) == dmer"TCG"d

    # Bounds checking
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 0)
    @test_throws BoundsError Base.setindex(dmer"TAGC"d, DNA_A, 5)
    @test_throws BoundsError Base.setindex(dmer""d, DNA_A, 1)

    # Test with larger bit integers
    d256 = DNAOligo{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_G, 3) == DNAOligo{UInt256}(dna"TAGCAT")
    @test Base.setindex(d256, DNA_C, 1) == DNAOligo{UInt256}(dna"CAGCAT")
end

@testset "capacity" begin
    # Create a zero-BPS alphabet for testing
    struct ZeroBPSAlphabet <: Alphabet end
    Base.eltype(::Type{ZeroBPSAlphabet}) = DNA
    BioSequences.encode(::ZeroBPSAlphabet, ::DNA) = zero(UInt)
    BioSequences.decode(::ZeroBPSAlphabet, ::Unsigned) = DNA_A
    BioSequences.BitsPerSymbol(::ZeroBPSAlphabet) = BioSequences.BitsPerSymbol{0}()

    # Test zero BPS: capacity should be clamped to typemax(Int)
    @test capacity(Oligo{ZeroBPSAlphabet, UInt8}) == clamp(typemax(UInt8), Int)
    @test capacity(Oligo{ZeroBPSAlphabet, UInt32}) == clamp(typemax(UInt32), Int)
    @test capacity(Oligo{ZeroBPSAlphabet, UInt128}) == typemax(Int)
    @test_throws ErrorException from_integer(
        Oligo{ZeroBPSAlphabet, UInt8}, zero(UInt8), 256
    )

    # Zero-BPS values have no coding bits, so from_integer must retain only the
    # requested length regardless of its integer input.
    zero_bps_type = Oligo{ZeroBPSAlphabet, UInt8}
    for value in [zero(UInt8), one(UInt8), typemax(UInt8)]
        @test from_integer(zero_bps_type, value, 0).x == 0x00
    end
    @test from_integer(zero_bps_type, typemax(UInt8), 1).x == 0x01
    zero_bps = from_integer(zero_bps_type, zero(UInt8), 4)

    symbols = fill(DNA_A, 4)
    zero_bps_kmer = Kmer{ZeroBPSAlphabet, 4}(symbols)
    @test length(zero_bps_kmer) == 4
    @test zero_bps_kmer.data === ()
    @test Oligo{ZeroBPSAlphabet, UInt8}(symbols) === zero_bps

    # Packed extraction must not read source storage for a zero-BPS alphabet.
    zero_bps_source = LongSequence{ZeroBPSAlphabet}(UInt64[], UInt(4))
    @test Kmers.extract_packed(
        Kmers.Copyable(),
        typeof(zero_bps_kmer),
        zero_bps_source,
        1,
        BioSequences.BitsPerSymbol{0}(),
    ) === zero_bps_kmer
    @test Kmers.extract_packed(
        Kmers.Copyable(),
        zero_bps_type,
        zero_bps_source,
        1,
        4,
        BioSequences.BitsPerSymbol{0}(),
    ) === zero_bps
    @test @inferred(
        Kmers.unsafe_extract(Kmers.GenericRecoding(), zero_bps_type, symbols, 1, 4)
    ) === zero_bps

    selected = @inferred zero_bps[[true, false, true, false]]
    @test typeof(selected) === zero_bps_type
    @test length(selected) == 2
    @test selected.x == 0x02

    # Test non-zero BPS: capacity should be in range 0:div(8 * sizeof(U), B)
    for (A, bps) in [(DNAAlphabet{2}, 2), (DNAAlphabet{4}, 4), (AminoAcidAlphabet, 8)]
        for U in [UInt8, UInt32, UInt64]
            cap = capacity(Oligo{A, U})
            max_possible = div(8 * sizeof(U), bps)
            @test 0 <= cap <= max_possible
        end
    end

    # Test that larger backing type gives larger or equal capacity
    @test capacity(DNAOligo{UInt32}) <= capacity(DNAOligo{UInt64})
    @test capacity(AAOligo{UInt32}) <= capacity(AAOligo{UInt64})
    @test capacity(DNAOligo{UInt64}) <= capacity(DNAOligo{UInt256})
    @test capacity(AAOligo{UInt128}) <= capacity(AAOligo{UInt512})

    # Instances forward to their concrete type's capacity.
    for oligomer in (
            DNAOligo{UInt8}("TAG"),
            Oligo{DNAAlphabet{4}, UInt32}("TAGN"),
            AAOligo{UInt64}("WPK"),
        )
        @test @inferred(capacity(oligomer)) == capacity(typeof(oligomer))
    end
end

@testset "translate" begin
    @testset "Output backing type" begin
        @test typeof(translate(DNAOligo{UInt8}(dna"ATG"))) === AAOligo{UInt16}
        @test typeof(translate(DNAOligo{UInt16}(dna"ATG"))) === AAOligo{UInt32}
        @test typeof(translate(DNAOligo{UInt32}(dna"ATG"))) === AAOligo{UInt64}
        @test typeof(translate(DNAOligo{UInt64}(dna"ATG"))) === AAOligo{UInt128}
        @test typeof(translate(DNAOligo{UInt128}(dna"ATG"))) === AAOligo{UInt256}
    end

    # Basic translation - 2-bit alphabets
    @testset "Basic 2-bit" begin
        # DNA 2-bit with different backing types
        @test string(translate(dmer"ATGTAA"d)) == "M*"
        @test string(translate(DNAOligo{UInt32}(dna"ATGTAA"))) == "M*"
        @test string(translate(DNAOligo{UInt64}(dna"ATGTAA"))) == "M*"

        # RNA 2-bit with different backing types
        @test string(translate(dmer"AUGUAA"r)) == "M*"
        @test string(translate(RNAOligo{UInt32}(rna"AUGUAA"))) == "M*"
        @test string(translate(RNAOligo{UInt64}(rna"AUGUAA"))) == "M*"

        # Longer sequences
        @test string(translate(dmer"TCTACACCCTAG"d)) == "STP*"
        @test string(translate(dmer"UCUACACCCUAG"r)) == "STP*"
        s = randdnaseq(249)
        t = translate(DNAOligo{UInt512}(s))
        @test LongSequence(t) == translate(s)
    end

    # Basic translation - 4-bit alphabets
    @testset "Basic 4-bit" begin
        # DNA 4-bit with different backing types
        d = Oligo{DNAAlphabet{4}, UInt64}(dna"TGGCCCGATTGA")
        @test string(translate(d)) == "WPD*"

        # UInt128 works with 4-bit since capacity is smaller
        d128 = Oligo{DNAAlphabet{4}, UInt128}(dna"ATGTAA")
        @test string(translate(d128)) == "M*"

        # RNA 4-bit with different backing types
        r = Oligo{RNAAlphabet{4}, UInt64}(rna"UGGCCCGAUUGA")
        @test string(translate(r)) == "WPD*"

        r32 = Oligo{RNAAlphabet{4}, UInt32}(rna"AUGUAA")
        @test string(translate(r32)) == "M*"

        r128 = Oligo{RNAAlphabet{4}, UInt128}(rna"AUGUAA")
        @test string(translate(r128)) == "M*"

        # Long sequence
        s = randdnaseq(252)
        t = translate(Oligo{DNAAlphabet{4}, UInt1024}(s))
        @test LongSequence(t) == translate(s)
    end

    # Different genetic codes
    @testset "Different genetic codes" begin
        # Vertebrate mitochondrial code: AGA and AGG are stop codons
        vert_mito = BioSequences.vertebrate_mitochondrial_genetic_code
        @test string(translate(dmer"ATGAGA"d; code = vert_mito)) == "M*"  # AGA is stop

        # Standard code: AGA and AGG code for R (Arginine)
        @test string(translate(dmer"ATGAGA"d)) == "MR"
        @test string(translate(dmer"ATGAGG"d)) == "MR"

        # Test with different backing types
        @test string(translate(DNAOligo{UInt64}(dna"ATGAGA"); code = vert_mito)) == "M*"
        @test string(translate(RNAOligo{UInt32}(rna"AUGAGA"); code = vert_mito)) == "M*"

        # Test with 4-bit alphabets
        @test string(translate(Oligo{DNAAlphabet{4}, UInt64}(dna"ATGAGA"); code = vert_mito)) == "M*"
    end

    # alternative_start flag
    @testset "alternative_start" begin
        # Without alternative_start: TTG codes for L (Leucine)
        @test string(translate(dmer"TTGCCC"d; alternative_start = false)) == "LP"
        @test string(translate(dmer"UUGCCC"r; alternative_start = false)) == "LP"

        # With alternative_start: first codon becomes M regardless
        @test string(translate(dmer"TTGCCC"d; alternative_start = true)) == "MP"
        @test string(translate(dmer"UUGCCC"r; alternative_start = true)) == "MP"
    end

    # Ambiguous codons (only for 4-bit alphabets)
    @testset "Ambiguous codons" begin
        # With allow_ambiguous_codons=true (default), ambiguous codons translate
        d_ambig = Oligo{DNAAlphabet{4}, UInt128}("AAAACWGCSWTARACADA")
        @test string(translate(d_ambig)) == "KTAJBX"

        # With allow_ambiguous_codons=false, ambiguous codons throw
        @test_throws Exception translate(d_ambig; allow_ambiguous_codons = false)

        # Test various ambiguous nucleotides
        # W = A or T, so TWG could be AAG (K) or TAG (*)
        # With allow_ambiguous, should give X (ambiguous)
        d_w = Oligo{DNAAlphabet{4}, UInt64}(dna"ATGTWG")
        result_w = translate(d_w; allow_ambiguous_codons = true)
        @test length(result_w) == 2  # M and something
    end

    # Error: Length not divisible by 3
    @testset "Length not divisible by 3" begin
        @test_throws ArgumentError translate(dmer"A"d)
        @test_throws ArgumentError translate(dmer"UG"r)
        @test_throws ArgumentError translate(dmer"CUGUAGUUGUCGC"r)
        @test_throws ArgumentError translate(Oligo{RNAAlphabet{4}, UInt32}(rna"AUGC"))
    end

    # Error: Sequences with gaps (only for 4-bit alphabets)
    @testset "Sequences with gaps" begin
        @test_throws Exception translate(Oligo{RNAAlphabet{4}, UInt64}(rna"-UGAUG"))
        @test_throws Exception translate(Oligo{DNAAlphabet{4}, UInt64}(dna"AT-ATG"))
        @test_throws Exception translate(Oligo{RNAAlphabet{4}, UInt64}(rna"AUGAU-"))
        @test_throws Exception translate(Oligo{DNAAlphabet{4}, UInt64}(dna"A--"))
    end

    # Error: Input type capacity too large for output
    @testset "Capacity overflow" begin
        # These tests have loaded BitIntegers, so the maximum backing integer
        # is currently UInt1024.
        # This means 1024-bit 2-bit alphabets overflow, but not 4-bit alphabets.
        @test_throws ErrorException translate(DNAOligo{UInt1024}(dna"ATG"))
        @test_throws ErrorException translate(RNAOligo{UInt1024}(rna"AUG"))

        # Even empty sequences should error due to type-based capacity check
        @test_throws ErrorException translate(DNAOligo{UInt1024}(dna""))
        @test_throws ErrorException translate(RNAOligo{UInt1024}(rna""))

        # Not an overflow - 512 bits fit
        @test string(translate(RNAOligo{UInt512}(rna""))) == ""

        # 4-bit nucleotides do not overflow even when 1024 bits
        @test string(translate(Oligo{RNAAlphabet{4}, UInt1024}(rna"AUG"))) == "M"
    end

    # Edge cases
    @testset "Edge cases" begin
        # Empty sequence (length 0 is divisible by 3, produces 0 AA)
        for A in (DNAAlphabet{2}, RNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{4}),
                U in (UInt8, UInt32), alternative_start in (false, true)
            input = Oligo{A, U}("")
            result = translate(input; alternative_start)
            empty_result = empty(typeof(result))
            @test length(result) == 0
            @test isempty(result)
            @test result == empty_result
            @test hash(result) == hash(empty_result)
            @test string(result) == ""
            @test result.x == zero(result.x)
        end

        # Very long sequence (near capacity) - 2-bit
        # For UInt64 with 2-bit DNA, capacity is 31, so 30 nt = 10 AA
        long_dna = dmer"ATGATGATGATGATGATGATGATGATG"d
        @test length(translate(long_dna)) == 9
        @test LongSequence(translate(long_dna)) == aa"MMMMMMMMM"
    end
end

@testset "Misc" begin
    d = AAOligo{UInt32}("WPK")
    @test only([d]') === d
end

@testset "Random extended-width Oligos" begin
    for T in (DNAOligo{UInt256}, Oligo{DNAAlphabet{4}, UInt512})
        oligomer = @inferred rand(StableRNG(SEED), T, capacity(T))
        @test typeof(oligomer) === T
        @test length(oligomer) == capacity(T)
    end
end
