module TestKmers

using Test
using Random
using StableRNGs
using Kmers
using BioSequences
using BioSymbols
using StringViews

const SEED = 0xccfb2d5055d8c990

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

struct GenericNucAlphabet <: NucleicAcidAlphabet{8} end
Base.eltype(::Type{GenericNucAlphabet}) = DNA
BioSequences.symbols(::GenericNucAlphabet) = symbols(DNAAlphabet{4}())
BioSequences.encode(::GenericNucAlphabet, c::DNA) = BioSequences.encode(DNAAlphabet{4}(), c)
BioSequences.decode(::GenericNucAlphabet, c::UInt) =
    BioSequences.decode(DNAAlphabet{4}(), c)
BioSequences.BitsPerSymbol(::GenericNucAlphabet) = BioSequences.BitsPerSymbol{8}()

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

    @testset "Bad parameters" begin
        @test_throws Exception Kmer{DNAAlphabet, :foo, 1}("C")
        @test_throws Exception Kmer{DNAAlphabet, -1, 0}("")
        @test_throws Exception Kmer{DNAAlphabet, 1, 2}("A")
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

        @test_throws Exception eval(:(mer"ATCGATAG"k)) # invalid flag
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
        @test mer""a == mer""a
        @test mer"KMNUPQCX"a == mer"KMNUPQCX"a
        @test mer"PKMNEA"a != mer"PKMNE"a
        @test mer"PKM"a != mer"PK"a
        @test mer"IUDHLDJVIPOEJKWE"a != mer"IUDHLDJVIPOEJKW"a
    end

    @testset "Ordering" begin
        @test mer"UGCAG"r > mer"CGCAG"r
        @test mer"TCGGAAG"d > mer"TCGGAAC"d
        @test mer"OEWPM"a > mer"OEWP"a
        @test mer"UGCGA"r > mer"TGAGA"d

        @test cmp(mer"TAGCTA"d, mer"TACCTA"d) == 1
        @test cmp(mer"TAC"d, mer"TAGCA"d) == -1
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
        @test_throws Exception mer"TAG"d == dna"TAG"
    end
end

@testset "As integer" begin
    u1 = as_integer(mer"TAGCGA"d)
    u2 = as_integer(mer"TAGCGC"d)
    @test u1 != u2
    @test u1 isa Unsigned
    @test u2 isa Unsigned

    d = Set{UInt8}()
    for s in Iterators.product(repeat(["ACGU"], 4)...)
        push!(d, as_integer(RNAKmer{4}(join(s))))
    end
    @test length(d) == 256

    u3 = as_integer(mer"TGATGCTGTAGTCGTGA"d)
    @test u3 isa Unsigned

    u4 = as_integer(mer"KWPLKWPHWLM"a)
    @test u4 isa Unsigned

    @test_throws ArgumentError as_integer(mer"AAAAAAAAAAAAAAAAAAAAAAAA"a)
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

    @testset "pop, pop_first" begin
        m = mer"LNPQ"a
        @test (m = pop(m)) == mer"LNP"a
        @test (m = pop(m)) == mer"LN"a
        @test (m = pop(m)) == mer"L"a
        @test (m = pop(m)) == AAKmer{0}("")
        @test_throws ArgumentError pop(m)

        @test pop(mer"MDFFIJFKL"a) === mer"MDFFIJFK"a

        m = mer"UAGC"r
        @test (m = pop_first(m)) == mer"AGC"r
        @test (m = pop_first(m)) == mer"GC"r
        @test (m = pop_first(m)) == mer"C"r
        @test (m = pop_first(m)) == mer""r
        @test pop_first(mer"PKWIKMPPAVYWA"a) == mer"KWIKMPPAVYWA"a
    end

    @testset "Setindex" begin
        mer = mer"PLQVAK"a
        setindex = Base.setindex
        @test setindex(mer, 3, AA_K) == mer"PLKVAK"a
        @test setindex(mer, 1, AA_R) == mer"RLQVAK"a
        @test setindex(mer, 6, AA_M) == mer"PLQVAM"a
        @test_throws BoundsError setindex(mer, 0, AA_K)
        @test_throws BoundsError setindex(mer, 7, AA_K)

        mer = mer"ATGTCGTGA"d
        @test setindex(mer, 1, DNA_T) == mer"TTGTCGTGA"d
        @test setindex(mer, 5, DNA_C) == mer"ATGTCGTGA"d
        @test setindex(mer, 5, DNA_A) == mer"ATGTAGTGA"d

        mer = mer"PLAKCVMARYKW"a
        @test setindex(mer, 10, AA_Q) == mer"PLAKCVMARQKW"a
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

        @test iscanonical(mer"AGCTAG"d)
        @test iscanonical(mer""d)
        @test iscanonical(mer"GCGAAC"d)
        @test iscanonical(mer"AATT"d)
        @test !iscanonical(mer"GGATGC"d)
        @test !iscanonical(mer"TCGTGA"d)
        @test !iscanonical(mer"TTGAA"d)
    end
end

@testset "Translation" begin
    @testset "Forward translation" begin
        # Empty
        @test translate(mer""r) == mer""a
        @test translate(mer""d) == mer""a
        @test translate(Kmer{DNAAlphabet{4}, 0}("")) == mer""a

        # Not divisible by 3
        @test_throws Exception translate(mer"U"r)
        @test_throws Exception translate(mer"UGCA"r)
        @test_throws Exception translate(mer"GUCGAUUGUC"r)

        # Containing gaps
        @test_throws Exception translate(Kmer{DNAAlphabet{4}, 6}("CTGA-C"))
        @test_throws Exception translate(Kmer{RNAAlphabet{4}, 3}("UC-"))

        # Invalid alphabet
        @test_throws Exception transate(mer"CCC"a)
        @test_throws Exception transate(Kmer{CharAlphabet, 3}("GGG"))

        # Compare to LongSequence
        for s in [
            rna"UCGUAGUUCGAUUCUAUGCUGUAGUGGCAA",
            rna"UCGUAGGCGUAUUGCGCAAAGCGC",
            rna"UGCUAGUGUUCGAAA",
            rna"UCGUUAGUAAAA",
        ]
            for A in [DNAAlphabet{4}, RNAAlphabet{2}, DNAAlphabet{2}, RNAAlphabet{4}]
                ss = LongSequence{A}(s)
                @test collect(translate(ss)) == collect(translate(Kmer{A, length(s)}(s)))
            end
        end

        for s in [
            rna"UGCUGAWKVUDUGWUGUDHUAGUGCNUBGKUGCMGGSWC",
            rna"UCGUAGUCKGUCGUYCUGAGGWUGCUGANNUGCUGA",
            rna"CAGGCCAGWGCUGSSSCUGSMGKYVUCUAS",
        ]
            for A in [DNAAlphabet{4}, RNAAlphabet{4}]
                ss = LongSequence{A}(s)
                @test collect(translate(ss)) == collect(translate(Kmer{A, length(s)}(s)))
            end
        end

        # Skip 1, the index of gap (which cannot be translated)
        A = alphabet(RNA)
        for i in 2:16, j in 2:16, k in 2:16
            mer = Kmer{RNAAlphabet{4}, 3}((A[i], A[j], A[k]))
            @test only(translate(mer)) == only(translate(LongSequence(mer)))
        end
    end

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

            s = isempty(cset) ? mer"AAA"r : first(cset)
            @test delete(cset, s) == delete!(copy(set), s)
            @test filter(i -> first(i) == DNA_A, cset) ==
                  filter(i -> first(i) == DNA_A, set)
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
            @test reverse_translate(aa) === codonset
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

        code = Kmers.rev_standard_genetic_code
        @test length(code) == 27
        @test collect(code) == map(0x00:UInt8(26)) do i
            aa = reinterpret(AminoAcid, i)
            aa => reverse_translate(aa)
        end
        @test_throws Exception code[AA_Gap]
    end

    @testset "Custom reverse genetic code" begin
        fw_code = BioSequences.pterobrachia_mitochondrial_genetic_code
        code = ReverseGeneticCode(fw_code)
        for (aa, set) in code
            for codon in set
                if aa ∈ (AA_O, AA_U, AA_B, AA_J, AA_X, AA_Z)
                    continue
                end
                @test only(translate(LongSequence(codon); code=fw_code)) === aa
            end
        end
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

            # Bad byte in ASCII
            s = "TAGTCGTAGPATGC"
            @test_throws BioSequences.EncodeError collect(FwDNAMers{3}(s))
        end

        # Unconvertible alphabet
        @testset "Unconvertible alphabet" begin
            @test_throws Exception iterate(FwKmers{DNAAlphabet{4}, 2}(aa"TAGTGCA"))
        end

        # GenericRecoding
        s = dna"TGATGTCGTAGTGAgtagtaCCA"
        it = FwKmers{GenericNucAlphabet, 8}(s)
        @test collect(it) ==
              [Kmer{GenericNucAlphabet, 8}(s[i:(i + 7)]) for i in 1:(length(s) - 7)]
    end

    @testset "FwRvIterator" begin
        function naive_fwrv(s::NucSeq, len::Integer)
            A = typeof(Alphabet(s))
            T = Kmers.derive_type(Kmer{A, len})
            [
                (T(s[i:(i + len - 1)]), T(reverse_complement(s[i:(i + len - 1)]))) for
                i in 1:(length(s) - len + 1)
            ]
        end

        for s in ["", "TGATGCTGTA", "TAT"]
            Ts = [DNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{2}, RNAAlphabet{4}]
            srcs = Any[LongSequence{T}(LongDNA{2}(s)) for T in Ts]
            dsts = copy(srcs)
            for dst in dsts
                srcs_ = push!(copy(srcs), String(typeof(dst)(LongDNA{2}(s))))
                for src in srcs_
                    @test collect(FwRvIterator{typeof(Alphabet(dst)), 4}(src)) ==
                          naive_fwrv(dst, 4)
                end
            end
        end
    end

    @testset "CanonicalKmers" begin
        @testset "Aliases" begin
            @test CanonicalKmers{DNAAlphabet{4}, 3}("TAGCTAGA") isa CanonicalKmers
            @test CanonicalDNAMers{8}("TAGCTAGA") isa CanonicalKmers
            @test CanonicalRNAMers{9}("TAGCTAGA") isa CanonicalKmers
        end

        @testset "Only nucleic acids" begin
            @test_throws Exception CanonicalKmers{AminoAcidAlphabet, 3}("UAGCTGA")
        end

        @testset "Iteration" begin
            for s in [
                dna"TAGCTAGGACA",
                rna"UAGUCGUGAGA",
                "TAGCTAGAGGA",
                collect(codeunits("ATGCGAGGA")),
            ]
                seq = LongDNA{2}(s)
                cns = [canonical(seq[i:(i + 4)]) for i in 1:(length(seq) - 4)]
                for A in [DNAAlphabet{2}, DNAAlphabet{4}]
                    it = CanonicalKmers{A, 5}(s)
                    @test collect(it) == [Kmer{A, 5, 1}(i) for i in cns]
                end
            end

            s = dna"TAGTCGTGATGATAGTCTGAATGTC"
            it = CanonicalKmers{GenericNucAlphabet, 6}(s)
            @test collect(it) == [
                canonical(Kmer{GenericNucAlphabet, 6}(s[i:(i + 5)])) for
                i in 1:(length(s) - 5)
            ]

            s = "TAGTGTCGATGATC"
            it1 = CanonicalKmers{DNAAlphabet{2}, 4, String}(s)
            it2 = CanonicalDNAMers{4}(s)
            @test collect(it1) == collect(it2)
        end
    end

    @testset "UnambiguousKmers" begin
        for s in [dna"TAGCWSAGACYWNACGCNACG--", rna"UAGUCYWUAGCNUAHAGC-GAUGAGC"]
            res = [(s[i:(i + 2)], i) for i in 1:(length(s) - 2)]
            filter!(i -> all(iscertain, first(i)), res)
            for A in [DNAAlphabet{2}, RNAAlphabet{2}]
                resA = [(Kmer{A, 3, 1}(i), j) for (i, j) in res]
                it = UnambiguousKmers{A, 3}(s)
                @test collect(it) == resA
            end

            A = s isa LongDNA ? DNAAlphabet{2} : RNAAlphabet{2}
            pos = [i:(i + 3) for i in 1:(lastindex(s) - 3)]
            filter!(pos) do rng
                all(iscertain, s[rng])
            end
            resA = [(Kmer{A, 4, 1}(s[i]), first(i)) for i in pos]
            it = UnambiguousKmers{A, 4}(string(s))
            @test collect(it) == resA
        end

        # Copyable
        s = LongDNA{2}("TATCGGATAGGCAAA")
        v = collect(UnambiguousRNAMers{4}(s))
        v2 = [(DNAKmer{4}(s[i:(i + 3)]), i) for i in 1:(length(s) - 3)]
        @test v == v2

        # GenericRecoding
        s = dna"TAGCTKAGAGGAGAACWSGCGAGA"
        it = UnambiguousKmers{DNAAlphabet{2}, 4}(s)
        v = collect(it)
        v2 = [
            (DNAKmer{4}(s[i:(i + 3)]), i) for
            i in 1:(length(s) - 3) if all(iscertain, s[i:(i + 3)])
        ]
        @test v == v2

        s = LongSequence{GenericNucAlphabet}(dna"TGATCGTAGATGwATGTC")
        it = UnambiguousKmers{DNAAlphabet{2}, 7, LongSequence{GenericNucAlphabet}}(s)
        it2 = UnambiguousDNAMers{7}(s)
        @test collect(it) == collect(it2)

        # Bad byte in ASCII
        s = "TAGTCGTAGPATGC"
        @test_throws BioSequences.EncodeError collect(UnambiguousDNAMers{3}(s))
    end

    @testset "SpacedKmers" begin
        function test_naive_spaced(A, seq, k, space)
            T = Kmers.derive_type(Kmer{A, k})
            v = [T(seq[i:(i + k - 1)]) for i in 1:space:(length(seq) - k + 1)]
            @test collect(SpacedKmers{A, k, space}(seq)) == v
        end

        for (s, A) in Any[
            ("TA-NGAKATCGAWTAGA", DNAAlphabet{4}),
            ("AUGCUGAUGAGUCGUAG", RNAAlphabet{2}),
            ("KLMYUPOKQMMNLVYRW", AminoAcidAlphabet),
        ]
            test_naive_spaced(A, s, 3, 2)
            test_naive_spaced(A, s, 2, 4)
            test_naive_spaced(A, codeunits(s), 3, 3)
        end
        test_naive_spaced(DNAAlphabet{2}, rna"UAGUCGUAGUAG", 4, 3)
        test_naive_spaced(RNAAlphabet{4}, dna"TAGCCWKMMNAGCTV", 2, 3)

        it = SpacedDNAMers{3, 4}("TAGAWWWW")
        @test_throws BioSequences.EncodeError collect(it)
    end

    @testset "Each codon" begin
        for s in Any[
            "TAGCGATAT",
            b"UAUGCUGAA",
            dna"TAGGCTATA",
            LongDNA{2}(dna"TAGCTAGAGGA"),
            rna"UGAUUCGUUGA",
            LongRNA{2}(rna"UAGUCGUGAGUA"),
        ]
            @test each_codon(DNA, s) === SpacedDNAMers{3, 3}(s)
            @test each_codon(RNA, s) === SpacedRNAMers{3, 3}(s)
            if s isa BioSequence{<:DNAAlphabet}
                @test each_codon(s) === SpacedDNAMers{3, 3}(s)
            elseif s isa BioSequence{<:RNAAlphabet}
                @test each_codon(s) === SpacedRNAMers{3, 3}(s)
            end
        end
    end
end

@testset "StringViews" begin
    s = StringView(collect(codeunits("ATGCTGATGATCGTATGATGTCGAAA")))
    it = FwRvIterator{DNAAlphabet{2}, 9}(s)
    @test collect(it) == map(1:(length(s) - 8)) do i
        a = DNAKmer{9}(s[i:(i + 8)])
        (a, reverse_complement(a))
    end
end

@testset "fx_hash" begin
    if UInt == UInt64
        for (kmer, h) in Any[
            (mer"TAG"a, 0x55dbbe22bb3e4a13),
            (mer"KPWAK"a, 0x10203d1c885b7467),
            (mer"TAGCTAG"d, 0xa76409341339d05a),
            (mer""a, 0x0000000000000000),
            (mer""r, 0x0000000000000000),
            (mer"UGAUGCA"r, 0xdd7c97ae4ca204b4),
        ]
            @test fx_hash(kmer) === h
        end
    end
end

@testset "Construction utils" begin
    @testset "Unsafe extract" begin
        seq = dna"TTGCTAGGGATTCGAGGATCCTCTAGAGCGCGGCACGATCTTAGCAC"
        unsafe_extract = Kmers.unsafe_extract
        @test unsafe_extract(Kmers.FourToTwo(), DNAKmer{6, 1}, seq, 3) ==
              DNAKmer{6}(seq[3:8])
        @test unsafe_extract(Kmers.FourToTwo(), DNAKmer{36, 2}, seq, 2) ==
              DNAKmer{36}(seq[2:37])

        seq = LongDNA{2}(seq)
        @test unsafe_extract(Kmers.TwoToFour(), Kmer{DNAAlphabet{4}, 6, 1}, seq, 3) ==
              Kmer{DNAAlphabet{4}, 6}(seq[3:8])
        @test unsafe_extract(Kmers.TwoToFour(), Kmer{DNAAlphabet{4}, 36, 3}, seq, 2) ==
              Kmer{DNAAlphabet{4}, 36}(seq[2:37])

        @test unsafe_extract(Kmers.Copyable(), DNAKmer{6, 1}, seq, 3) ==
              DNAKmer{6}(seq[3:8])
        @test unsafe_extract(Kmers.Copyable(), DNAKmer{36, 2}, seq, 2) ==
              DNAKmer{36}(seq[2:37])

        seq = codeunits(String(seq))
        @test unsafe_extract(Kmers.AsciiEncode(), DNAKmer{6, 1}, seq, 3) ==
              DNAKmer{6}(seq[3:8])
        @test unsafe_extract(Kmers.AsciiEncode(), DNAKmer{36, 2}, seq, 2) ==
              DNAKmer{36}(seq[2:37])

        seq = LongSequence{CharAlphabet}("中国¨Å!人大æ网")
        @test unsafe_extract(Kmers.GenericRecoding(), Kmer{CharAlphabet, 3, 2}, seq, 4) ==
              Kmer{CharAlphabet, 3}("Å!人")
    end

    @testset "Unsafe shift from" begin
        ushift = Kmers.unsafe_shift_from

        seq = dna"TTGCTAGGGATTCGAGGATCCTCTAGAGCGCGGCACGATCTTAGCAC"
        mer = Kmer{DNAAlphabet{4}, 9}("TAGwKwADH")
        @test ushift(Kmers.Copyable(), mer, seq, 4, Val(3)) ==
              Kmer{DNAAlphabet{4}, 9}("wKwADHCTA")

        mer = mer"TAGCATCG"d
        @test ushift(Kmers.FourToTwo(), mer, seq, 4, Val(3)) == mer"CATCGCTA"d

        seq = LongDNA{2}(seq)
        mer = Kmer{DNAAlphabet{4}, 9}("TAGwKwADH")
        @test ushift(Kmers.TwoToFour(), mer, seq, 2, Val(3)) ==
              Kmer{DNAAlphabet{4}, 9}("wKwADHTGC")

        seq = codeunits(String(seq))
        mer = mer"KWPLCVAKVM"a
        @test ushift(Kmers.AsciiEncode(), mer, seq, 5, Val(4)) == mer"CVAKVMTAGG"a

        seq = LongSequence{CharAlphabet}("中国¨Å!人大æ网")
        mer = Kmer{CharAlphabet, 5, 3}("中国¨Å!")
        @test ushift(Kmers.GenericRecoding(), mer, seq, 6, Val(3)) ==
              Kmer{CharAlphabet, 5, 3}("Å!人大æ")
    end
end

@testset "Random kmers" begin
    @testset "Complete alphabets" begin
        @test rand(StableRNG(SEED), DNAKmer{10}) == mer"GATAAACTTG"d
        @test rand(StableRNG(SEED), DNAKmer{10, 1}) == mer"GATAAACTTG"d
        @test rand(StableRNG(SEED), DNAKmer{0}) == mer""d
        @test rand(StableRNG(SEED), DNAKmer{41}) == mer"TTACGTTCAGGGGCAGTCGAGATCGGCTCCGGATAAACTTG"d
        @test rand(StableRNG(SEED), RNAKmer{4}) == mer"CUUG"r
        @test rand(StableRNG(SEED), DNAKmer{4}) == mer"CTTG"d
        @test rand(StableRNG(SEED), DNAKmer{64}) == mer"GGGGCAGTCGAGATCGGCTCCGGATAAACTTGCATACAAGGCGTAGAGAAGAGTTTTACGTTCA"d
    end

    @testset "Incomplete alphabets" begin
        @test rand(StableRNG(SEED), AAKmer{10}) == mer"WFTEYAFNSW"a
        @test rand(StableRNG(SEED), AAKmer{0, 0}) == mer""a
        @test rand(StableRNG(SEED), AAKmer{16}) == mer"WFTEYAFNSWVMQHLS"a
    end

    @testset "Four-bit alphabets" begin
        @test rand(StableRNG(SEED), Kmer{DNAAlphabet{4}, 12}) == Kmer{DNAAlphabet{4}, 12}("GGCGCTGAAATG")
        @test rand(StableRNG(SEED), Kmer{RNAAlphabet{4}, 12}) == Kmer{RNAAlphabet{4}, 12}("GGCGCUGAAAUG")
        @test rand(StableRNG(SEED), Kmer{RNAAlphabet{4}, 33}) == Kmer{RNAAlphabet{4}, 33}("GGGACGGCGCUGAAAUGAAAGCCGGAACUAGUA")
        @test rand(StableRNG(SEED), Kmer{RNAAlphabet{4}, 32}) == Kmer{RNAAlphabet{4}, 32}("AAAGCCGGAACUAGUAGGACGGCGCUGAAAUG")
    end

    @testset "Custom alphabets" begin
        @test rand(StableRNG(SEED), Kmer{CharAlphabet, 7}) == Kmer{CharAlphabet, 7}("~t2\$0\x037")
        @test rand(StableRNG(SEED), Kmer{CharAlphabet, 0}) == Kmer{CharAlphabet, 0}("")

        @test rand(StableRNG(SEED), Kmer{CharAlphabet, 51}) == Kmer{CharAlphabet, 51}(
            "~t2\$0\x0373@\$K8/Ryuz\x144\x18x\e\x11\nmW&[{\tO7XY'O\r?\0c2P\n=\x03^)Bu\x02s"
        )
    end

    @testset "Instances" begin
        @test rand(StableRNG(SEED), mer"PKWSJMVTYWB"a) == AA_J
    end
end

end # module
