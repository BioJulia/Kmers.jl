@testset "Translation" begin

sampler = BioSequences.SamplerWeighted(
    dna"ACGTMRSVWYHKDBN",
    vcat(fill(0.225, 4), fill(0.00909, 10))
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
                    translate(seq, alternative_start=alternative) ==
                    translate(kmer, alternative_start=alternative)
                )          
            end
        end
    end
end

# Throws when ambiguous
@test_throws Exception translate(
    Kmer{RNAAlphabet{4}}("AUGCCGCMA"),
    allow_ambiguous_codons=false
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
