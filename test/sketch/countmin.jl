@testset "CountMinSketch tests" begin
    @testset "CountMinSketch normal constructor" begin
        cms = CountMinSketch{UInt16}(4, 100)
        @test eltype(cms.sketch) == UInt16
        @test size(cms.sketch) == (4, 100)

        # table too small @test_throws Exception CountMinSketch{UInt8}(4,0)
        # too few tables
        @test_throws Exception CountMinSketch{UInt8}(0,100)
        @test_throws Exception CountMinSketch{UInt8}(100,100)
    end

    @testset "CountMinSketch size/eltype" begin
        cms = CountMinSketch{UInt16}(4, 100)
        @test eltype(cms.sketch) == UInt16
        @test eltype(cms) == UInt16
        @test size(cms.sketch) == (4, 100)
        @test size(cms) == (4, 100)
    end

    @testset "CountMinSketch push/pop/add" begin
        cms = CountMinSketch{UInt16}(4, 100)
        @test cms[hash(1)] == UInt16(0)
        push!(cms, hash(1))
        @test cms[hash(1)] == UInt16(1)
        push!(cms, hash(1))
        @test cms[hash(1)] == UInt16(2)
        pop!(cms, hash(1))
        @test cms[hash(1)] == UInt16(1)
        add!(cms, hash(1), 100)
        @test cms[hash(1)] == UInt16(101)
        add!(cms, hash(1), -10)
        @test cms[hash(1)] == UInt16(91)
    end

    @testset "CountMinSketch getindex" begin
        cms = CountMinSketch{UInt16}(4, 100)
        item = "ABCD"
        @test cms[item] == UInt16(0)
        push!(cms, item)
        @test cms[item] == UInt16(1)
        push!(cms, item)
        @test cms[item] == UInt16(2)
        pop!(cms, item)
        @test cms[item] == UInt16(1)
        add!(cms, item, 100)
        @test cms[item] == UInt16(101)
        add!(cms, item, -10)
        @test cms[item] == UInt16(91)
    end

    @testset "CountMinSketch over/underflow" begin
        cms = CountMinSketch{UInt8}(4, 100)
        item = "ABCD"

        @test cms[item] == 0
        add!(cms, item, 256)
        @test cms[item] == 255
        add!(cms, item, -256)
        @test cms[item] == 0
    end

    @testset "collisionrate(CountMinSketch)" begin
        cms = CountMinSketch{UInt16}(4, 100)

        # Fill the first half of each table. This is super dodgy
        for table in 1:4
            for i in 1:50
                cms.sketch[table, i] = 1
            end
        end
        @test collisionrate(cms) == 0.5^4
    end
end # testset CountMinSketch

