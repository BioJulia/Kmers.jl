@testset "Translation" begin
    sampler = BioSequences.SamplerWeighted(
        dna"ACGTMRSVWYHKDBN",
        vcat(fill(0.225, 4), fill(0.00909, 10)),
    )

    for A in (RNAAlphabet, DNAAlphabet)
        for N in (2, 4)
            for len in [3, 15, 33, 66]
                for alternative in (true, false)
                    seq = if N == 2
                        randseq(A{2}(), len)
                    else
                        randseq(A{4}(), sampler, len)
                    end
                    kmer = Kmer{A{N}}(seq)
                    @test (
                        translate(seq; alternative_start = alternative) ==
                            translate(kmer; alternative_start = alternative)
                    )
                end
            end
        end
    end

    # Throws when ambiguous
    @test_throws Exception translate(
        Kmer{RNAAlphabet{4}}("AUGCCGCMA"),
        allow_ambiguous_codons = false,
    )

    # Not divisible by 3
    @test_throws Exception translate(mer"UG"r)
    @test_throws Exception translate(mer"TAGCTTAA"d)
    @test_throws Exception translate(mer"CUGUAGUUGUCGC"r)
    @test_throws Exception translate(mer"AGCGA"d)

    # Cannot transla AA seq
    @test_throws MethodError translate(mer"LLVM"aa)
    @test_throws MethodError translate(mer"ATG"aa)
end # translation

@testset "CodonSet" begin
    CodonSet = Kmers.CodonSet

    SAMPLE_SOURCES = Any[
        [mer"UAG"r, mer"ACC"r, mer"ACC"r, mer"UGG"r],
        RNACodon[],
        [mer"AAA"r, mer"ACC"r, mer"AAA"r, mer"UCA"r, mer"UCC"r],
        (i for i in (mer"AGC"r, mer"AGA"r, mer"UUU"r)),
        (mer"AAC"r, mer"AGG"r),
        (mer"UUG"r,),
    ]

    @testset "Construction and basics" begin
        @test isempty(CodonSet())

        # Constuct the sets and basic properties
        for codons in SAMPLE_SOURCES
            set = Set(codons)
            codonset = CodonSet(codons)
            @test issetequal(set, codonset)
            @test length(codonset) == length(set)
        end

        # Fails with non-codons
        @test_throws MethodError CodonSet([(RNA_A, RNA_G)])
        @test_throws MethodError CodonSet((mer"UA"r,))
        @test_throws MethodError CodonSet([rna"AGG", rna"GGG"])
        @test_throws MethodError CodonSet([1, 2, 3])
    end

    SAMPLE_CODONSETS = map(CodonSet, SAMPLE_SOURCES)

    @testset "Iteration" begin
        for things in SAMPLE_SOURCES
            @test sort!(collect(CodonSet(things))) == sort!(collect(Set(things)))
        end

        @test iterate(CodonSet()) === nothing
        codonset = CodonSet((mer"UUU"r,))
        codon, state = iterate(codonset)
        @test codon == mer"UUU"r
        @test iterate(codonset, state) === nothing
    end

    @testset "Membership" begin
        codonset = CodonSet([mer"ACC"r, mer"UAG"r, mer"UUU"r])
        @test mer"ACC"r in codonset
        @test mer"UAG"r in codonset
        @test mer"UUU"r in codonset
        @test !in(mer"GAA"r, codonset)
        @test !in(mer"AAA"r, codonset)
    end

    @testset "Modifying" begin
        # Push
        s1 = CodonSet([mer"GGA"r, mer"UGU"r])
        s2 = push(s1, mer"GGA"r)
        @test s1 == s2
        s3 = push(s2, mer"GAG"r)
        @test Set(s3) == Set([mer"GGA"r, mer"UGU"r, mer"GGA"r, mer"GAG"r])

        # Delete
        s4 = delete(s3, mer"GAG"r)
        @test s2 == s4
        s5 = delete(s4, mer"UGU"r)
        @test only(s5) == mer"GGA"r
        s6 = delete(s5, mer"UUU"r)
        @test s5 == s6
        s7 = delete(s6, mer"GGA"r)
        @test isempty(s7)
    end

    @testset "Set operations" begin
        for c1 in SAMPLE_CODONSETS, c2 in SAMPLE_CODONSETS
            s1, s2 = Set(c1), Set(c2)
            for operation in [union, intersect, setdiff, symdiff]
                @test Set(operation(c1, c2)) == operation(s1, s2)
            end
            @test issubset(c1, c2) == issubset(s1, s2)
        end
    end

    @testset "Filter" begin
        predicates = [
            (i -> i[2] == RNA_G),
            (i -> isodd(length(i))), # always true for codons
            (i -> i[1] == i[3]),
            (i -> i[2] != RNA_A),
        ]
        for codonset in SAMPLE_CODONSETS, predicate in predicates
            @test Set(filter(predicate, codonset)) == filter(predicate, Set(codonset))
        end
    end
end # CodonSet

@testset "Reverse translation" begin
    CodonSet = Kmers.CodonSet
    code2 = ReverseGeneticCode(BioSequences.trematode_mitochondrial_genetic_code)
    for (rvcode, fwcode) in [
            (Kmers.rev_standard_genetic_code, BioSequences.standard_genetic_code),
            (code2, BioSequences.trematode_mitochondrial_genetic_code),
        ]
        @test reverse_translate(aa"", rvcode) == CodonSet[]
        observed = Dict{AminoAcid, CodonSet}()
        for aa in symbols(AminoAcidAlphabet())
            # Gap cannot be reverse translated
            aa in (AA_Gap,) && continue
            observed[aa] = only(reverse_translate(LongAA([aa]), rvcode))
        end

        # Length and iteration of ReverseGeneticCode
        @test length(rvcode) == length(symbols(AminoAcidAlphabet())) - 1 # all but AA_Gap
        @test sort!(collect(rvcode); by = first) == sort!(collect(observed); by = first)

        flipped = Dict(v => k for (k, v) in observed)
        for (codonset, aa) in flipped
            # Ambigous AA are reverse translated to a set, each of which directly
            # translate to a specific non-ambig AA, so this is not reversible
            BioSequences.isambiguous(aa) && continue
            # AA_O and AA_U are specifically in the reverse standard genetic code
            # but not in the genetic code.
            # This is because they can be unambiguously revtranslated, but not
            # forward translated.
            aa in (AA_O, AA_U) && continue
            for codon in codonset
                @test only(translate(LongRNA{2}(codon); code = fwcode)) == aa
            end
        end

        # Special AA
        @test only(reverse_translate(aa"O", rvcode)) == CodonSet((mer"UAG"r,))
        @test only(reverse_translate(aa"U", rvcode)) == CodonSet((mer"UGA"r,))

        # Ambiguous amino acids
        for (ambig, elements) in [
                (AA_J, [AA_I, AA_L]),
                (AA_Z, [AA_E, AA_Q]),
                (AA_B, [AA_D, AA_N]),
                (
                    AA_X,
                    [
                        AA_A,
                        AA_R,
                        AA_N,
                        AA_D,
                        AA_C,
                        AA_Q,
                        AA_E,
                        AA_G,
                        AA_H,
                        AA_I,
                        AA_L,
                        AA_K,
                        AA_M,
                        AA_F,
                        AA_P,
                        AA_S,
                        AA_T,
                        AA_W,
                        AA_Y,
                        AA_V,
                    ],
                ),
            ]
            c1 = only(reverse_translate(LongAA([ambig]), rvcode))
            c2 = foldl(elements; init = CodonSet()) do old, aa
                union(old, reverse_translate(aa, rvcode))
            end
            @test c1 == c2
        end

        # Test error on gap
        @test_throws Exception reverse_translate(aa"-")
    end
end # reverse translation
