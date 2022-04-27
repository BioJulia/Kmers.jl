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

@testset "SpacedKmers" begin
    @testset "SpacedKmers DNA" begin
        s = randdnaseq(500)
        s2 = LongDNA{2}(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) == collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) == length(SpacedKmers(s2, Val{31}(), 50)) == 10
        
        @test collect(SpacedKmers(s, Val{201}(), 50)) == collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) == length(SpacedKmers(s2, Val{201}(), 50)) == 6
    end
    
    @testset "SpacedKmers RNA" begin
        s = randrnaseq(500)
        s2 = LongRNA{2}(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) == collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) == length(SpacedKmers(s2, Val{31}(), 50)) == 10
        
        @test collect(SpacedKmers(s, Val{201}(), 50)) == collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) == length(SpacedKmers(s2, Val{201}(), 50)) == 6
    end
    
    @testset "SpacedKmers AA" begin
        s = randaaseq(500)
        s2 = LongAA(s)
        @test collect(SpacedKmers(s, Val{31}(), 50)) == collect(SpacedKmers(s2, Val{31}(), 50))
        @test length(SpacedKmers(s, Val{31}(), 50)) == length(SpacedKmers(s2, Val{31}(), 50)) == 10
    
        @test collect(SpacedKmers(s, Val{201}(), 50)) == collect(SpacedKmers(s2, Val{201}(), 50))
        @test length(SpacedKmers(s, Val{201}(), 50)) == length(SpacedKmers(s2, Val{201}(), 50)) == 6
    end
end