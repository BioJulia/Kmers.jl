@testset "De Bruijn Neighbors" begin
    @test collect(fw_neighbors(DNAKmer("ACG")))  == map(DNAKmer, ["CGA",  "CGC",  "CGG",  "CGT"])
    @test collect(fw_neighbors(DNAKmer("GGGG"))) == map(DNAKmer, ["GGGA", "GGGC", "GGGG", "GGGT"])
    @test collect(fw_neighbors(RNAKmer("ACG")))  == map(RNAKmer, ["CGA",  "CGC",  "CGG",  "CGU"])
    @test collect(fw_neighbors(RNAKmer("GGGG"))) == map(RNAKmer, ["GGGA", "GGGC", "GGGG", "GGGU"])
end
