@testset "Access and Iterations" begin
    dna_kmer = mer"ACTG"dna
    rna_kmer = mer"ACUG"rna
    aa_kmer  = mer"MVXN"aa

    @testset "Access DNA Kmer" begin
        @test dna_kmer[1] == DNA_A
        @test dna_kmer[2] == DNA_C
        @test dna_kmer[3] == DNA_T
        @test dna_kmer[4] == DNA_G
        
        @test dna_kmer[1:3] == mer"ACT"dna
        @test dna_kmer[2:4] == mer"CTG"dna

        # Access indexes out of bounds
        @test_throws BoundsError dna_kmer[-1]
        @test_throws BoundsError dna_kmer[0]
        @test_throws BoundsError dna_kmer[5]
        @test_throws BoundsError getindex(dna_kmer,-1)
        @test_throws BoundsError getindex(dna_kmer, 0)
        @test_throws BoundsError getindex(dna_kmer, 5)
        @test_throws BoundsError dna_kmer[3:7]
    end

    @testset "Iteration through DNA Kmer" begin
        @test iterate(DNAKmer("ACTG")) == (DNA_A, 2)

        @test iterate(DNAKmer("ACTG"), 1) == (DNA_A, 2)
        @test iterate(DNAKmer("ACTG"), 4) == (DNA_G, 5)

        @test iterate(DNAKmer("ACTG"), 1)  !== nothing
        @test iterate(DNAKmer("ACTG"), 4)  !== nothing
        @test iterate(DNAKmer("ACTG"), 5)  === nothing
        @test_throws BoundsError iterate(DNAKmer("ACTG"), -1)

        dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
        @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_kmer)])
    end

    @testset "Access RNA Kmer" begin
        @test rna_kmer[1] == RNA_A
        @test rna_kmer[2] == RNA_C
        @test rna_kmer[3] == RNA_U
        @test rna_kmer[4] == RNA_G
        
        @test rna_kmer[1:3] == mer"ACU"rna
        @test rna_kmer[2:4] == mer"CUG"rna

        # Access indexes out of bounds
        @test_throws BoundsError rna_kmer[-1]
        @test_throws BoundsError rna_kmer[0]
        @test_throws BoundsError rna_kmer[5]
        @test_throws BoundsError getindex(rna_kmer, -1)
        @test_throws BoundsError getindex(rna_kmer, 0)
        @test_throws BoundsError getindex(rna_kmer, 5)
        @test_throws BoundsError rna_kmer[3:7]
    end

    @testset "Iteration through RNA Kmer" begin
        @test iterate(RNAKmer("ACUG")) == (RNA_A, 2)
        

        @test iterate(RNAKmer("ACUG"), 1) == (RNA_A, 2)
        @test iterate(RNAKmer("ACUG"), 4) == (RNA_G, 5)
        

        @test iterate(RNAKmer("ACUG"), 1)  !== nothing
        @test iterate(RNAKmer("ACUG"), 4)  !== nothing
        @test iterate(RNAKmer("ACUG"), 5)  === nothing
        @test_throws BoundsError iterate(RNAKmer("ACUG"), -1)

        rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
        @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_kmer)])
    end
    
    @testset "Access AA Kmer" begin
        @test aa_kmer[1] == AA_M
        @test aa_kmer[2] == AA_V
        @test aa_kmer[3] == AA_X
        @test aa_kmer[4] == AA_N

        @test aa_kmer[1:3] == mer"MVX"aa
        @test aa_kmer[2:4] == mer"VXN"aa

        # Access indexes out of bounds
        @test_throws BoundsError aa_kmer[-1]
        @test_throws BoundsError aa_kmer[0]
        @test_throws BoundsError aa_kmer[5]
        @test_throws BoundsError getindex(aa_kmer,-1)
        @test_throws BoundsError getindex(aa_kmer, 0)
        @test_throws BoundsError getindex(aa_kmer, 5)
        @test_throws BoundsError aa_kmer[3:7]
    end
end
