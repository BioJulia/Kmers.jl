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
    @testset "Scalar indexing" begin
        m = mer"TGATGCTAGTAGTATTCTATAG"d
        @test m isa Mer{22}

        @test m[1] == first(m) == DNA_T
        @test m[3] == DNA_A
        @test last(m) == m[22] == DNA_G

        @test_throws BoundsError m[0]
        @test_throws BoundsError m[-1]
        @test_throws BoundsError m[23]

        for s in
            [dna"TAGCAAC", dna"TWKKSVVDNA-A", rna"UGUGUCA", rna"UGUCGWS", aa"PLLKMDDSH"]
            m = Kmer{typeof(Alphabet(s)), length(s)}(s)
            @test first(m) == first(s)
            @test last(m) == last(s)
            for i in [1, 3, 5]
                @test s[i] == m[i]
            end
        end

        # Weirdly, this throws ArgumentError
        @test_throws Exception first(DNAKmer{0}(""))
        @test_throws Exception last(RNAKmer{0}(""))
    end

    @testset "Unit ranges" begin
        m = mer"POKDGTWDIKVL"a
        @test m isa Mer{12}

        @test m[1:3] === mer"POK"a
        @test m[2:6] === mer"OKDGT"a
        @test m[6:(end - 1)] === mer"TWDIKV"a

        @test m[eachindex(m)] === m
        @test m[Base.OneTo(4)] === mer"POKD"a

        @test_throws BoundsError m[0:4]
        @test m[0:-1] == AAKmer{0}("")
        @test_throws BoundsError m[2:13]
    end

    @testset "With vector of indices" begin
        m = mer"UGCUGAUCGUAU"r
        @test m isa Mer{12}

        @test m[[1, 3, 5]] == mer"UCG"r
        @test m[[12, 9, 7]] == mer"UGU"r
        @test m[Int[]] == Kmer{RNAAlphabet{2}, 0}("")

        @test_throws BoundsError m[[2, 8, 15]]
        @test_throws BoundsError m[[0, 1]]
        @test_throws BoundsError m[[13]]
    end

    @testset "Logical indexing" begin
        m = Kmer{CharAlphabet, 4}("ØÆGD")
        @test m[[true, false, true, false]] == Kmer{CharAlphabet, 2}("ØG")
        @test m[trues(4)] === m
        @test m[falses(4)] === Kmer{CharAlphabet, 0}("")

        @test_throws BoundsError m[[true, false, true, true, true]]
        @test_throws BoundsError m[[false, false, true]]
        @test_throws BoundsError m[trues(5)]
    end
end

@testset "Modification" begin
    @testset "push, push_first" begin
        m = mer"UHALSAP"a
        @test push(m, AA_W) == mer"UHALSAPW"a
        @test push(push(m, AA_W), AA_M) === mer"UHALSAPWM"a
        @test push_first(m, AA_Gap) == mer"-UHALSAP"a
        @test push_first(push(m, AA_K), AA_H) == mer"HUHALSAPK"a

        @test push(m, 'K') == mer"UHALSAPK"a
        @test push(mer"TAG"d, RNA_A) == mer"TAGA"d
    end

    @testset "shift, shiftfirst" begin
        m = mer"PDOFPOLEF"a
        v = collect(m)
        for aa in aa"PLLMWFVB"
            m = shift(m, aa)
            @test m isa Mer{9}
            popfirst!(push!(v, aa))
            @test collect(m) == v
        end

        m = mer"AUGCGUA"r
        v = collect(m)
        for dna in dna"TAGTGTGCTA"
            m = shift_first(m, dna)
            @test m isa Mer{7}
            pop!(pushfirst!(v, dna))
            @test collect(m) == v
        end
    end

    @testset "pop" begin
        m = mer"LNPQ"a
        @test (m = pop(m)) == mer"LNP"a
        @test (m = pop(m)) == mer"LN"a
        @test (m = pop(m)) == mer"L"a
        @test (m = pop(m)) == AAKmer{0}("")
        @test_throws ArgumentError pop(m)

        @test pop(mer"MDFFIJFKL"a) === mer"MDFFIJFK"a
    end
end

@testset "Biological operations" begin
    for s in [
        dna"",
        aa"",
        LongDNA{2}(dna"TAGTGCA"),
        LongRNA{2}(rna"UGCUGUAA"),
        dna"TGASWKHVAAN--A",
        rna"UAGUCUYMNS",
        aa"LKHWSYYVQN",
        LongSequence{CharAlphabet}("LKDSJ"),
        LongSequence{CharAlphabet}("κ𝚶⊸∑Γ"),
    ]
        m = Kmer{typeof(Alphabet(s)), length(s)}(s)

        # Reverse
        @test collect(reverse(m)) == reverse(collect(m))
        @test collect(reverse(m)) == collect(reverse(s))

        # The rest of the operations are only for nucleotides
        isa(Alphabet(s), NucleicAcidAlphabet) || continue

        # Complement
        @test collect(complement(s)) == collect(complement(m))

        # Reverse complement
        rv = reverse_complement(m)
        @test collect(reverse_complement(s)) == collect(rv)

        # Canonical
        can = canonical(m)
        @test collect(can) == collect(canonical(s))
        @test can ≤ m
        if can === m
            @test m ≤ rv
        else
            @test can === rv
            @test rv ≤ m
        end
    end
end

@testset "Translation" begin
    @testset "CodonSet" begin
        codons = [RNACodon((i, j, k)) for i in mer"UACG"r, j in mer"UACG"r, k in mer"UACG"r]
        @test length(Set(codons)) == 64
        sources = [[], codons[[1, 4, 8]], codons, codons[rand(Bool, 64)], codons[[4, 8]]]
        csets = map(CodonSet, sources)
        sets = map(Set, sources)

        # Basic properties
        for (cset, set) in zip(csets, sets)
            @test cset == set
            @test sort!(collect(cset)) == sort!(collect(set))
            @test length(cset) == length(set)
            @test isempty(cset) == isempty(set)
            for i in set
                @test i ∈ cset
            end
        end

        for (si, ci) in zip(sets, csets), (sj, cj) in zip(sets, csets)
            @test issubset(si, sj) == issubset(ci, cj)
            for f in [union, setdiff, intersect, symdiff]
                @test Set(f(ci, cj)) == f(si, sj)
            end
        end
    end

    @testset "Standard reverse genetic code" begin
        seq = LongAA(collect((i for i in alphabet(AminoAcid) if i ∉ (AA_Gap, AA_U, AA_O))))
        codonsets = reverse_translate(seq)
        seen_codons = Set{RNACodon}()
        for (codonset, aa) in zip(codonsets, seq)
            if isambiguous(aa)
                bits = zero(compatbits(aa))
                for codon in codonset
                    bits |= compatbits(only(translate(codon)))
                end
                # selenocysteine and Pyrrolysine have bits
                # 0x00300000. However, translating normal
                # codons cannot get these amino acids,
                # so we ignore them by masking their bits
                @test bits == (compatbits(aa) & 0x000fffff)
            else
                @test isdisjoint(seen_codons, codonset)
                union!(seen_codons, codonset)
                for codon in codonset
                    @test only(translate(codon)) == aa
                end
            end
        end
        @test length(seen_codons) == 64
    end

    @testset "Custom reverse genetic code" begin
        # TODO!
    end
end

@testset "Printing" begin
    function test_print(s, str)
        @test string(s) == str
        io = IOBuffer()
        print(io, s)
        @test String(take!(io)) == str
    end

    for s in [
        dna"",
        aa"",
        LongDNA{2}(dna"TAGTGCA"),
        LongRNA{2}(rna"UGCUGUAA"),
        dna"TGASWKHVAAN--A",
        rna"UAGUCUYMNS",
        aa"LKHWSYYVQN",
    ]
        test_print(s, string(s))
    end
end

@testset "Iterators" begin
    @testset "Forward iteration" begin
        @testset "Aliases" begin
            @test FwKmers{DNAAlphabet{2}, 3}(dna"TAGA") isa
                  FwKmers{DNAAlphabet{2}, 3, LongDNA{4}}
            @test FwDNAMers{4}(rna"UAGC") isa FwKmers{DNAAlphabet{2}, 4, LongRNA{4}}
            @test FwRNAMers{4}(dna"TACA") isa FwKmers{RNAAlphabet{2}, 4, LongDNA{4}}
            @test FwAAMers{4}(aa"LKCY") isa FwKmers{AminoAcidAlphabet, 4, LongAA}
        end

        @testset "Smaller than K" begin
            @test isempty(FwDNAMers{3}(dna"TA"))
            @test isempty(FwAAMers{9}(aa"AOPJVPES"))
            @test isempty(FwKmers{RNAAlphabet{4}, 6}(dna"ATGGA"))
        end

        @testset "Conversible alphabets" begin
            for (seqs, alphabets) in [
                (
                    [LongDNA{2}("TGATGGCGTAGTA"), LongRNA{2}("UCGUGCUA"), LongDNA{2}("")],
                    [DNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{2}, RNAAlphabet{4}],
                ), # From two-bit
                (
                    [dna"TAGTCTGAC", rna"UAGUCGAUUAGGCC"],
                    [DNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{2}, RNAAlphabet{4}],
                ), # From four-bit
            ]
                for seq in seqs, alphabet in alphabets
                    v1 = collect(FwKmers{alphabet, 3}(seq))
                    v2 = [Kmer{alphabet, 3, 1}(seq[i:(i + 2)]) for i in 1:(length(seq) - 2)]
                    @test v1 == v2
                end
            end
            for seq in [dna"TGWSNVNTGA", rna"C-GGAU-WSNUCG"]
                @test_throws Exception first(FwDNAMers{3}(seq))
                @test_throws Exception first(FwRNAMers{3}(seq))
            end
        end

        @testset "Four to two bit" begin
            for seq in [dna"TATGCTTCGTAGTCGTCGTTGCTA"]
                for seqq in [seq, LongRNA{4}(seq)]
                    filtered = typeof(seqq)([i for i in seqq if !isambiguous(i)])
                    for A in [DNAAlphabet{2}, RNAAlphabet{2}]
                        v1 = collect(FwKmers{A, 4}(seqq))
                        v2 = [
                            Kmer{A, 4, 1}(filtered[i:(i + 3)]) for
                            i in 1:(length(filtered) - 3)
                        ]
                        @test v1 == v2
                    end
                end
            end
        end

        @testset "From ASCII bytes" begin
            str = "TaghWS-TGnADbkWWMSTV"
            T = FwKmers{DNAAlphabet{4}, 4}
            mers = collect(T(str))
            for source in
                [str, view(str, 1:lastindex(str)), codeunits(str), Vector(codeunits(str))]
                @test collect(T(source)) == mers
            end
        end

        # Unconvertible alphabet
        # TODO
        # Error in the constructor?
        #@test_throws FwDNAMers{}
    end
end

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
