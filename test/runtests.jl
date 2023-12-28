module TestKmers

using Test
using Random
using Kmers
using BioSequences
using BioSymbols

include("utils.jl")

@testset "BioSequences Interface" begin
    for A in
        [DNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{2}, RNAAlphabet{4}, AminoAcidAlphabet]
        for K in (1, 9, 116)
            @test BioSequences.has_interface(
                BioSequence,
                Kmers.derive_type(Kmer{A, K}),
                rand(collect(A()), K),
                false,
            )
        end
    end
end

struct CharSymbol <: BioSymbol
    x::Char
end
BioSymbols.prefix(::CharSymbol) = "Char"
BioSymbols.type_text(::CharSymbol) = "CharSymbol"
BioSymbols.isterm(::CharSymbol) = false

# These two are not interface, but useful for the tests below
# (e.g. converting the sequence to a string)
Base.convert(::Type{Char}, x::CharSymbol) = x.x
Base.Char(x::CharSymbol) = convert(Char, x)
Base.convert(::Type{CharSymbol}, x::Char) = CharSymbol(x)
BioSymbols.isgap(::CharSymbol) = false

# TODO: Should BioSymbols be updated to remove this?
BioSymbols.isvalid(::CharSymbol) = true

struct CharAlphabet <: Alphabet end
Base.eltype(::Type{CharAlphabet}) = CharSymbol
BioSequences.symbols(::CharAlphabet) = ntuple(i -> CharSymbol(Char(i - 1)), Val{128}())
BioSequences.encode(::CharAlphabet, c::CharSymbol) = reinterpret(UInt32, c.x) % UInt
BioSequences.decode(::CharAlphabet, c::UInt) = CharSymbol(reinterpret(Char, c % UInt32))
BioSequences.BitsPerSymbol(::CharAlphabet) = BioSequences.BitsPerSymbol{32}()

const ALPHABETS = [
    DNAAlphabet{2}(),
    RNAAlphabet{2}(),
    DNAAlphabet{4}(),
    RNAAlphabet{4}(),
    AminoAcidAlphabet(),
    CharAlphabet(),
]

@testset "Construction" begin
    # Fundamentals
    dna = dna"TAGCTAAC"
    mer = Kmer{DNAAlphabet{2}, length(dna)}(dna)
    @test mer isa Kmer{DNAAlphabet{2}, length(dna)}
    @test DNACodon == DNAKmer{3, 1}
    @test RNACodon == RNAKmer{3, 1}

    for A in ALPHABETS
        Ta = typeof(A)
        for L in [0, 3, 11, 41]
            for i in 1:3
                # Fundamentals and length
                s = randseq(A, SamplerUniform(symbols(A)), L)
                mer = Kmer{Ta, L}(collect(s))
                @test mer isa Kmer{Ta, L}
                @test length(mer) == L
                @test string(mer) == string(s)

                # From string
                mer2 = Kmer{Ta, L}(string(s))
                @test mer == mer2
                @test mer === mer2
                @test string(mer) == string(mer2)

                # From LongSequence
                mer3 = Kmer{Ta, L}(s)
                @test mer === mer3
            end
        end
    end

    # Construct from string
    @testset "Construct from string" begin
        for s in ["TAG", "ACCGAGCC", "TGATGCTATTAGG"]
            L = length(s)
            for ss in (s, view(s, 1:lastindex(s)))
                @test DNAKmer{L, 1}(ss) == DNAKmer{L}(ss)
                @test string(DNAKmer{L}(ss)) == ss
            end
        end

        s = "UHVKALRIQURPFLSMOF"
        @test string(AAKmer{18}(s)) == s

        for s in ["αβγδϵ", "", "中国人大网"]
            L = length(s)
            for ss in (s, view(s, 1:lastindex(s)))
                sq = Kmer{CharAlphabet, L}(ss)
                @test string(sq) == s
                @test [Char(i) for i in sq] == collect(ss)
            end
        end

        # Wrong length - also for iterators of unknown size
    end

    @testset "Wrong length" begin
        @test_throws Exception DNAKmer{4}("TAC")
        @test_throws Exception DNAKmer{4}("TACCA")
        @test_throws Exception Kmer{CharAlphabet, 2}(['T'])
        @test_throws Exception AAMer{3}((AminoAcid(i) for i in "WPLK" if true))
    end

    @testset "Length must be given explicitly" begin
        for s in ["TACA", ""]
            @test_throws Exception DNAKmer("TACGA")
            @test string(DNAKmer{length(s)}(s)) == s
        end
        @test_throws Exception AAMer(aa"WPLKM")
        @test collect(AAKmer{5}(aa"WPLKM")) == collect(aa"WPLKM")
    end

    @testset "Kmer literal" begin
        @test collect(mer"TGAGTCA"d) == collect(dna"TGAGTCA")
        @test collect(mer"WQOPMKAP"a) == collect(aa"WQOPMKAP")
        @test collect(mer"UAUCGGAUC"r) == collect(rna"UAUCGGAUC")
    end

    @testset "Construct from Biosequences" begin
        @testset "Construct from LongSequence" begin
            for seq in [
                dna"TAGGCA",
                rna"UUCUGUGAGUCC",
                aa"TTCGGAA",
                LongSequence{CharAlphabet}("HELLO"),
            ]
                for sq in [seq, view(seq, 2:lastindex(seq))]
                    A = typeof(Alphabet(sq))
                    @test Kmer{A, length(sq)}(sq) == Kmer{A, length(sq)}(string(sq))
                    @test string(Kmer{A, length(sq)}(sq)) == string(sq)
                    @test_throws Exception Kmer{A, length(sq) + 1}(sq)
                end
            end
        end

        @testset "Construct from kmer" begin
            m = mer"TAGCGTTA"d
            m2 = DNAKmer{8}(m)
            @test m === m2
            @test_throws Exception DNAKmer{7}(m)
            m3 = RNAKmer{8}(m)
            @test m3 === mer"UAGCGUUA"r
            @test_throws Exception RNAKmer{9}(m)
            @test_throws Exception AAKmer{8}(m)
        end

        # From generic biosequence - TODO
    end

    @testset "Construct from iterable" begin
        m1 = DNAKmer{6}((i for i in dna"GCGATC"))
        m2 = DNAKmer{6}((i for i in dna"ATCGATGCAA" if i ∈ (DNA_A, DNA_C)))
        @test m1 === mer"GCGATC"d
        @test m2 === mer"ACACAA"d
        m3 = DNAKmer{4}((i for i in rna"GAUC" if true))
        @test m3 === mer"GATC"d
    end
end

@testset "Comparison" begin
    @testset "Equality" begin
        @test mer"KMNUPQCX"a == mer"KMNUPQCX"a
        @test mer"PKMNEA"a != mer"PKMNE"a
        @test mer"IUDHLDJVIPOEJKWE"a != mer"IUDHLDJVIPOEJKW"a
    end

    @testset "Ordering" begin
        @test mer"UGCAG"r > mer"CGCAG"r
        @test mer"TCGGAAG"d > mer"TCGGAAC"d
        @test mer"OEWPM"a > mer"OEWP"a
        @test mer"UGCGA"r > mer"TGAGA"d
    end

    @testset "Hashing, isless and isequal" begin
        @test hash(mer"POSMDGF"a, UInt(15)) === hash(mer"POSMDGF"a, UInt(15))
        @test isequal(mer"POSMDGF"a, mer"POSMDGF"a)

        # Same, but DNA/RNA
        m1 = mer"TAGCTA"d
        m2 = mer"UAGCUA"r
        @test isequal(m1, m2)
        @test hash(m1) === hash(m2)
        m3 = Kmer{DNAAlphabet{4}}(m1)
        m4 = Kmer{RNAAlphabet{4}}(m2)
        @test isequal(m3, m4)
        @test hash(m3) === hash(m4)

        # Other kmers
        # This throws because we want kmer hashing to be maximally fast,
        # which implies they must have a different hashing strategy from
        # other BioSequences, which implies they can't be isequal
        @test_throws Exception isequal(mer"UGCUGA"r, mer"UGCUGA"a)
        @test !isequal(mer"UGCAC"r, mer"UGCGA"r)

        # Other sequences
        @test_throws Exception dna"TAG" == mer"TAG"d
    end
end

@testset "Access" begin
    # Scalar indexing

    # Index with UnitRange

    # Index with vector of indices

    # Boolean indexing
end

@testset "Modification" begin
    # Push, pushfirst

    # Shift, shift_first

    # Pop
end

@testset "Biological operations" begin
    # Reverse

    # Complement

    # Reverse complement

    # Canonical
end

@testset "Translation" begin end

@testset "Iterators" begin end

# include("construction_and_conversion.jl")
# include("comparisons.jl")
# include("length.jl")
# include("access.jl")
# include("random.jl")
# include("find.jl")
# include("print.jl")
# include("transformations.jl")
# include("mismatches.jl")
# include("debruijn_neighbors.jl")
# include("iteration.jl")
# include("translation.jl")
#include("shuffle.jl")

end # module
