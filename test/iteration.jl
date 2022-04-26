@testset "EveryKmer" begin
    @testset "EveryKmer DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)
        @test collect(EveryKmer(s, Val{31}())) == collect(EveryKmer(s2, Val{31}()))
        @test length(EveryKmer(s, Val{31}())) == length(EveryKmer(s2, Val{31}())) == 470
        
        @test collect(EveryKmer(s, Val{201}())) == collect(EveryKmer(s2, Val{201}()))
        @test length(EveryKmer(s, Val{201}())) == length(EveryKmer(s2, Val{201}())) == 300
    end
    
    @testset "EveryKmer RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)
        @test collect(EveryKmer(s, Val{31}())) == collect(EveryKmer(s2, Val{31}()))
        @test length(EveryKmer(s, Val{31}())) == length(EveryKmer(s2, Val{31}())) == 470
        
        @test collect(EveryKmer(s, Val{201}())) == collect(EveryKmer(s2, Val{201}()))
        @test length(EveryKmer(s, Val{201}())) == length(EveryKmer(s2, Val{201}())) == 300
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