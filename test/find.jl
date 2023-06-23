@testset "Find" begin
    kmer = DNAKmer("ACGAG")

    @test findnext(DNA_A, kmer, 1) == 1
    @test findnext(DNA_C, kmer, 1) == 2
    @test findnext(DNA_G, kmer, 1) == 3
    @test findnext(DNA_T, kmer, 1) == nothing
    @test findnext(DNA_A, kmer, 2) == 4

    @test_throws BoundsError findnext(DNA_A, kmer, 0)
    @test findnext(DNA_A, kmer, 6) === nothing

    @test findprev(DNA_A, kmer, 5) == 4
    @test findprev(DNA_C, kmer, 5) == 2
    @test findprev(DNA_G, kmer, 5) == 5
    @test findprev(DNA_T, kmer, 5) == nothing
    @test findprev(DNA_G, kmer, 4) == 3

    @test findprev(DNA_A, kmer, 0) === nothing
    @test_throws BoundsError findprev(DNA_A, kmer, 6)

    @test findfirst(DNA_A, kmer) == 1
    @test findfirst(DNA_G, kmer) == 3
    @test findlast(DNA_A, kmer) == 4
    @test findlast(DNA_G, kmer) == 5

    kmer = AAKmer("AMVKFPSMT")

    @test findnext(AA_A, kmer, 1) == 1
    @test findnext(AA_M, kmer, 1) == 2
    @test findnext(AA_V, kmer, 1) == 3
    @test findnext(AA_K, kmer, 1) == 4
    @test findnext(AA_F, kmer, 1) == 5
    @test findnext(AA_P, kmer, 1) == 6
    @test findnext(AA_S, kmer, 1) == 7
    @test findnext(AA_M, kmer, 1) == 2
    @test findnext(AA_T, kmer, 1) == 9

    @test findnext(AA_F, kmer, 4) == 5
    @test findprev(AA_F, kmer, 4) == nothing
    @test findnext(AA_A, kmer, 7) == nothing
    @test findnext(AA_M, kmer, 5) == 8

    @test findfirst(AA_M, kmer) == 2
    @test findlast(AA_M, kmer) == 8
end
