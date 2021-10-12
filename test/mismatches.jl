@testset "Mismatches" begin
    function test_mismatches(a, b)
        count = 0
        for (x, y) in zip(a, b)
            count += x != y
        end
        @test mismatches(a, b) === mismatches(b, a) === count
    end

    for len in 1:64, _ in 1:10
        a = random_dna_kmer(len)
        b = random_dna_kmer(len)
        test_mismatches(DNAKmer(a), DNAKmer(b))
        test_mismatches(Kmer{DNAAlphabet{4}}(a), Kmer{DNAAlphabet{4}}(b))
        
        a = random_rna_kmer(len)
        b = random_rna_kmer(len)
        test_mismatches(RNAKmer(a), RNAKmer(b))
        test_mismatches(Kmer{RNAAlphabet{4}}(a), Kmer{RNAAlphabet{4}}(b))
        
        a = AAKmer(random_aa(len))
        b = AAKmer(random_aa(len))
        test_mismatches(a, b)
    end
end

@testset "Matches" begin
    function test_matches(a, b)
        count = 0
        for (x, y) in zip(a, b)
            count += x == y
        end
        @test matches(a, b) === matches(b, a) === count
    end

    for len in 1:64, _ in 1:10
        a = random_dna_kmer(len)
        b = random_dna_kmer(len)
        test_matches(DNAKmer(a), DNAKmer(b))
        test_matches(Kmer{DNAAlphabet{4}}(a), Kmer{DNAAlphabet{4}}(b))
        
        a = random_rna_kmer(len)
        b = random_rna_kmer(len)
        test_matches(RNAKmer(a), RNAKmer(b))
        test_matches(Kmer{RNAAlphabet{4}}(a), Kmer{RNAAlphabet{4}}(b))
        
        a = AAKmer(random_aa(len))
        b = AAKmer(random_aa(len))
        test_matches(a, b)
    end
end