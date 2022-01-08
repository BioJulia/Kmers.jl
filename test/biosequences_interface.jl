@testset "BioSequences Interface" begin
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{DNAAlphabet{2},31}), rand(ACGT, 31), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{DNAAlphabet{4},31}), rand(ACGT, 31), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{RNAAlphabet{2},31}), rand(ACGU, 31), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{RNAAlphabet{4},31}), rand(ACGU, 31), false)
    
    
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{DNAAlphabet{2},200}), rand(ACGT, 200), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{DNAAlphabet{4},200}), rand(ACGT, 200), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{RNAAlphabet{2},200}), rand(ACGU, 200), false)
    @test BioSequences.has_interface(BioSequence, Kmers.kmertype(Kmer{RNAAlphabet{4},200}), rand(ACGU, 200), false)
end