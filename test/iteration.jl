@testset "EveryKmer" begin
    @testset "EveryKmer DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)
        # Kmer and sequence Alphabets match.
        @test collect(EveryKmer(s, Val{31}())) == collect(EveryKmer(s2, Val{31}()))
        @test length(EveryKmer(s, Val{31}())) == length(EveryKmer(s2, Val{31}())) == 470

        @test collect(EveryKmer(s, Val{201}())) == collect(EveryKmer(s2, Val{201}()))
        @test length(EveryKmer(s, Val{201}())) == length(EveryKmer(s2, Val{201}())) == 300

        # Kmer and sequence Alphabets mismatch.
        s3 = dna"AC-TGAG--TGC"
        @test collect(EveryKmer{DNACodon}(s3)) == [
            (UInt64(4), Kmer(DNA_T, DNA_G, DNA_A)),
            (UInt64(5), Kmer(DNA_G, DNA_A, DNA_G)),
            (UInt64(10), Kmer(DNA_T, DNA_G, DNA_C)),
        ]
    end

    @testset "EveryKmer RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)
        # Kmer and sequence Alphabets match.
        @test collect(EveryKmer(s, Val{31}())) == collect(EveryKmer(s2, Val{31}()))
        @test length(EveryKmer(s, Val{31}())) == length(EveryKmer(s2, Val{31}())) == 470

        @test collect(EveryKmer(s, Val{201}())) == collect(EveryKmer(s2, Val{201}()))
        @test length(EveryKmer(s, Val{201}())) == length(EveryKmer(s2, Val{201}())) == 300

        # Kmer and sequence Alphabets mismatch.
        s3 = rna"AC-UGAG--UGC"
        @test collect(EveryKmer{RNACodon}(s3)) == [
            (UInt64(4), Kmer(RNA_U, RNA_G, RNA_A)),
            (UInt64(5), Kmer(RNA_G, RNA_A, RNA_G)),
            (UInt64(10), Kmer(RNA_U, RNA_G, RNA_C)),
        ]
    end

    @testset "EveryKmer AA" begin
        s = randaaseq(500)
        s2 = LongAA(s)
        @test collect(EveryKmer(s, Val{31}())) == collect(EveryKmer(s2, Val{31}()))
        @test length(EveryKmer(s, Val{31}())) == length(EveryKmer(s2, Val{31}())) == 470

        @test collect(EveryKmer(s, Val{201}())) == collect(EveryKmer(s2, Val{201}()))
        @test length(EveryKmer(s, Val{201}())) == length(EveryKmer(s2, Val{201}())) == 300
    end
end

@testset "SpacedKmers" begin
    @testset "SpacedKmers DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) ==
              collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) ==
              length(SpacedKmers(s2, Val{31}(), 50)) ==
              10

        @test collect(SpacedKmers(s, Val{201}(), 50)) ==
              collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) ==
              length(SpacedKmers(s2, Val{201}(), 50)) ==
              6

        s3 = dna"AC-TGAG--TGC"
        @test collect(SpacedKmers{DNACodon}(s3, 3)) == [
            (UInt64(4), Kmer(DNA_T, DNA_G, DNA_A)),
            (UInt64(10), Kmer(DNA_T, DNA_G, DNA_C)),
        ]
    end

    @testset "SpacedKmers RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) ==
              collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) ==
              length(SpacedKmers(s2, Val{31}(), 50)) ==
              10

        @test collect(SpacedKmers(s, Val{201}(), 50)) ==
              collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) ==
              length(SpacedKmers(s2, Val{201}(), 50)) ==
              6

        s3 = rna"AC-UGAG--UGC"
        @test collect(SpacedKmers{RNACodon}(s3, 3)) == [
            (UInt64(4), Kmer(RNA_U, RNA_G, RNA_A)),
            (UInt64(10), Kmer(RNA_U, RNA_G, RNA_C)),
        ]
    end

    @testset "SpacedKmers AA" begin
        s = randaaseq(500)
        s2 = LongAA(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) ==
              collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) ==
              length(SpacedKmers(s2, Val{31}(), 50)) ==
              10

        @test collect(SpacedKmers(s, Val{201}(), 50)) ==
              collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) ==
              length(SpacedKmers(s2, Val{201}(), 50)) ==
              6
    end
end

@testset "EveryCanonicalKmer" begin
    @testset "EveryCanonicalKmer DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)

        # Iterator generates expected results...
        ## 2-Bit DNA
        @test [(x[1], canonical(x[2])) for x in EveryKmer(s2, Val{31}())] == collect(EveryCanonicalKmer(s2, Val{31}()))

        @test [(x[1], canonical(x[2])) for x in EveryKmer(s2, Val{201}())] == collect(EveryCanonicalKmer(s2, Val{201}()))

        ## 4-Bit DNA
        @test [(x[1], canonical(x[2])) for x in EveryKmer(s, Val{31}())] == collect(EveryCanonicalKmer(s, Val{31}()))

        @test [(x[1], canonical(x[2])) for x in EveryKmer(s, Val{201}())] == collect(EveryCanonicalKmer(s, Val{201}()))

        # Test equivalency between different levels of bit compression...
        @test [x[2] for x in EveryCanonicalKmer(s, Val{31}())] == [x[2] for x in EveryCanonicalKmer(s2, Val{31}())]
        @test all(iscanonical.([x[2] for x in EveryCanonicalKmer(s, Val{31}())])) && all(iscanonical.([x[2] for x in EveryCanonicalKmer(s2, Val{31}())]))

        @test [x[2] for x in EveryCanonicalKmer(s, Val{201}())] == [x[2] for x in EveryCanonicalKmer(s2, Val{201}())]
        @test all(iscanonical.([x[2] for x in EveryCanonicalKmer(s, Val{201}())])) && all(iscanonical.([x[2] for x in EveryCanonicalKmer(s2, Val{201}())]))

        # Kmer and sequence Alphabets mismatch.
        s3 = dna"AC-TGAG--TGC"
        @test collect(EveryCanonicalKmer{DNACodon}(s3)) == [
            (UInt64(4), canonical(Kmer(DNA_T, DNA_G, DNA_A))),
            (UInt64(5), canonical(Kmer(DNA_G, DNA_A, DNA_G))),
            (UInt64(10), canonical(Kmer(DNA_T, DNA_G, DNA_C))),
        ]
    end

    @testset "EveryCanonicalKmer RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)

        # Iterator generates expected results...
        ## 2-Bit DNA
        @test [(x[1], canonical(x[2])) for x in EveryKmer(s2, Val{31}())] == collect(EveryCanonicalKmer(s2, Val{31}()))

        @test [(x[1], canonical(x[2])) for x in EveryKmer(s2, Val{201}())] == collect(EveryCanonicalKmer(s2, Val{201}()))

        ## 4-Bit DNA
        @test [(x[1], canonical(x[2])) for x in EveryKmer(s, Val{31}())] == collect(EveryCanonicalKmer(s, Val{31}()))

        @test [(x[1], canonical(x[2])) for x in EveryKmer(s, Val{201}())] == collect(EveryCanonicalKmer(s, Val{201}()))

        # Test equivalency between different levels of bit compression...
        @test [x[2] for x in EveryCanonicalKmer(s, Val{31}())] == [x[2] for x in EveryCanonicalKmer(s2, Val{31}())]
        @test all(iscanonical.([x[2] for x in EveryCanonicalKmer(s, Val{31}())])) && all(iscanonical.([x[2] for x in EveryCanonicalKmer(s2, Val{31}())]))

        @test [x[2] for x in EveryCanonicalKmer(s, Val{201}())] == [x[2] for x in EveryCanonicalKmer(s2, Val{201}())]
        @test all(iscanonical.([x[2] for x in EveryCanonicalKmer(s, Val{201}())])) && all(iscanonical.([x[2] for x in EveryCanonicalKmer(s2, Val{201}())]))

        s3 = rna"AC-UGAG--UGC"
        @test collect(EveryCanonicalKmer{RNACodon}(s3)) == [
            (UInt64(4), canonical(Kmer(RNA_U, RNA_G, RNA_A))),
            (UInt64(5), canonical(Kmer(RNA_G, RNA_A, RNA_G))),
            (UInt64(10), canonical(Kmer(RNA_U, RNA_G, RNA_C))),
        ]
    end
end

@testset "SpacedCanonicalKmers" begin
    @testset "SpacedCanonicalKmers DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s, Val{31}(), 50))
        @test collect(SpacedCanonicalKmers(s, Val{31}(), 50)) ==
              collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test length(SpacedCanonicalKmers(s, Val{31}(), 50)) ==
              length(SpacedCanonicalKmers(s2, Val{31}(), 50)) ==
              10

        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s, Val{201}(), 50))
        @test collect(SpacedCanonicalKmers(s, Val{201}(), 50)) ==
              collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test length(SpacedCanonicalKmers(s, Val{201}(), 50)) ==
              length(SpacedCanonicalKmers(s2, Val{201}(), 50)) ==
              6

        s3 = dna"AC-TGAG--TGC"
        @test collect(SpacedCanonicalKmers{DNACodon}(s3, 3)) == [
            (UInt64(4), canonical(Kmer(DNA_T, DNA_C, DNA_A))),
            (UInt64(10), canonical(Kmer(DNA_T, DNA_G, DNA_C))),
        ]
    end

    @testset "SpacedCanonicalKmers RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{31}(), 50)] == collect(SpacedCanonicalKmers(s, Val{31}(), 50))
        @test collect(SpacedCanonicalKmers(s, Val{31}(), 50)) ==
              collect(SpacedCanonicalKmers(s2, Val{31}(), 50))
        @test length(SpacedCanonicalKmers(s, Val{31}(), 50)) ==
              length(SpacedCanonicalKmers(s2, Val{31}(), 50)) ==
              10

        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test [(x[1], canonical(x[2])) for x in SpacedKmers(s2, Val{201}(), 50)] == collect(SpacedCanonicalKmers(s, Val{201}(), 50))
        @test collect(SpacedCanonicalKmers(s, Val{201}(), 50)) ==
              collect(SpacedCanonicalKmers(s2, Val{201}(), 50))
        @test length(SpacedCanonicalKmers(s, Val{201}(), 50)) ==
              length(SpacedCanonicalKmers(s2, Val{201}(), 50)) ==
              6

        s3 = rna"AC-UGAG--UGC"
        @test collect(SpacedCanonicalKmers{RNACodon}(s3, 3)) == [
            (UInt64(4), canonical(Kmer(RNA_U, RNA_C, RNA_A))),
            (UInt64(10), canonical(Kmer(RNA_U, RNA_G, RNA_C))),
        ]
    end
end
