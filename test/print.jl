@testset "Print" begin
    buf = IOBuffer()

    print(buf, DNAKmer("ACGT"))
    @test String(take!(buf)) == "ACGT"

    print(buf, RNAKmer("ACGU"))
    @test String(take!(buf)) == "ACGU"

    print(buf, Kmer{DNAAlphabet{4}}("ACGT"))
    @test String(take!(buf)) == "ACGT"

    print(buf, Kmer{RNAAlphabet{4}}("ACGU"))
    @test String(take!(buf)) == "ACGU"

    print(buf, AAKmer("AMVKFPSMT"))
    @test String(take!(buf)) == "AMVKFPSMT"
end

@testset "Show" begin
    buf = IOBuffer()

    show(buf, DNAKmer("AGAGT"))
    @test String(take!(buf)) == "AGAGT"

    show(buf, RNAKmer("AGAGU"))
    @test String(take!(buf)) == "AGAGU"

    show(buf, Kmer{DNAAlphabet{4}}("AGAGT"))
    @test String(take!(buf)) == "AGAGT"

    show(buf, Kmer{RNAAlphabet{4}}("AGAGU"))
    @test String(take!(buf)) == "AGAGU"

    print(buf, AAKmer("AMVKFPSMT"))
    @test String(take!(buf)) == "AMVKFPSMT"
end
