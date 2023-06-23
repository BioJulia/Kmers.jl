@testset "Transformations" begin
    function test_reverse(T, seq)
        revseq = reverse(T(seq))
        @test String(revseq) == reverse(seq)
    end

    function test_dna_complement(T, seq)
        comp = complement(T(seq))
        @test String(comp) == dna_complement(seq)
    end

    function test_rna_complement(T, seq)
        comp = complement(T(seq))
        @test String(comp) == rna_complement(seq)
    end

    function test_dna_revcomp(T, seq)
        revcomp = reverse_complement(T(seq))
        @test String(revcomp) == reverse(dna_complement(seq))
    end

    function test_rna_revcomp(T, seq)
        revcomp = reverse_complement(T(seq))
        @test String(revcomp) == reverse(rna_complement(seq))
    end

    @testset "Reverse" begin
        for len in 1:64, _ in 1:10
            test_reverse(DNAKmer{len}, random_dna_kmer(len))
            test_reverse(RNAKmer{len}, random_rna_kmer(len))
        end

        seq = dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAGATAC"
        @test reverse(seq[(length(seq) - 9):length(seq)]) == dna"CATAGAAAAA"
    end

    @testset "Complement" begin
        for len in 1:64, _ in 1:10
            test_dna_complement(DNAKmer{len}, random_dna_kmer(len))
            test_rna_complement(RNAKmer{len}, random_rna_kmer(len))
        end
    end

    @testset "Reverse Complement" begin
        for len in 1:64, _ in 1:10
            test_dna_revcomp(DNAKmer{len}, random_dna_kmer(len))
            test_rna_revcomp(RNAKmer{len}, random_rna_kmer(len))
        end
    end

    @testset "Canonical" begin
        @test canonical(DNAKmer{4, 1}("ACCG")) == DNAKmer{4, 1}("ACCG")
        @test canonical(DNAKmer{4, 1}("GCAC")) == DNAKmer{4, 1}("GCAC")
        @test canonical(RNAKmer{4, 1}("AAUU")) == RNAKmer{4, 1}("AAUU")
        @test canonical(RNAKmer{4, 1}("UGGA")) == RNAKmer{4, 1}("UCCA")
        @test canonical(RNAKmer{4, 1}("CGAU")) == RNAKmer{4, 1}("AUCG")
        @test canonical(RNAKmer{4, 1}("UGGA")) == RNAKmer{4, 1}("UCCA")
        @test canonical(DNAKmer{4, 1}("GCAC")) == DNAKmer{4, 1}("GCAC")
    end
end
