global reps = 10

@testset "Construction and Conversions" begin
    @test Kmer(DNA_A, DNA_G, DNA_T) === Kmer("AGT")
    @test Kmer(RNA_A, RNA_G, RNA_U) === Kmer("AGU")
    #@test Kmer(AA_R, AA_D, AA_C, AA_B) === Kmer("RDCB")

    @test DNAKmer(DNA_G, DNA_C, DNA_T) == Kmer("GCT")
    @test RNAKmer(RNA_G, RNA_U, RNA_C, RNA_U) == Kmer("GUCU")
    
    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(::Type{T}, seq::AbstractString) where {T<:Kmer}
        return String(T(seq)) == uppercase(seq)
    end
    
    # Check that RNAKmers can be constructed from a LongRNASeq
    #   LongSequence{A} → Kmer{A,K,N} → LongSequence{A}
    function check_longsequence_construction(::Type{T}, seq::S) where {T<:Kmer,S<:LongSequence}
        return S(T(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(::Type{T}, seq::LongSequence) where {T<:Kmer}
        return LongSequence(T(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{T} → Kmer → Vector{T}
    function check_nucarray_kmer(::Type{M}, seq::Vector{T}) where {T,M<:Kmer}
        return String([convert(Char, c) for c in seq]) == String(M(seq))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(::Type{T}, A2, seq::AbstractString) where {T<:Kmer}
        return String(LongSequence{A2}(T(LongSequence{A2}(seq)))) == uppercase(seq)
    end
    
    #=
    function check_uint_conversion(::Type{T}) where {T<:Kmer}
        U = BioSequences.encoded_data_type(T)
        uint = rand(typemin(U):U(one(U) << 2BioSequences.ksize(T) - 1))
        return convert(U, T(uint)) === uint
    end
    =#

    @testset "Kmer conversion" begin
        for len in [1, 16, 32, 64, 128]
            # String construction
            #   Check that kmers in strings survive round trip conversion:
            #   String → Kmer → String
            @test all(Bool[check_string_construction(DNAKmer{len},             random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Kmer{DNAAlphabet{4},len}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(RNAKmer{len},             random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Kmer{RNAAlphabet{4},len}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(AAKmer{len},              random_aa(len))       for _ in 1:reps])
            
            # Long(DNA|RNA)Seq Constructions
            #   Check that DNAKmers can be constructed from a Long(DNA|RNA)Seq
            #   Long(DNA|RNA)Seq → Kmer → Long(DNA|RNA)Seq
            @test all(Bool[check_longsequence_construction(Kmer{DNAAlphabet{2},len}, LongDNA{2}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{DNAAlphabet{4},len}, LongDNA{4}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{DNAAlphabet{4},len}, LongDNA{2}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{DNAAlphabet{2},len}, LongDNA{4}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{2},len}, LongRNA{2}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{4},len}, LongRNA{4}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{4},len}, LongRNA{2}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{2},len}, LongRNA{4}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(AAKmer{len},              LongAA(random_aa(len)))    for _ in 1:reps])
            
            # Check Kmer{A1}(::BioSequence{A2}) for compatible A1 and A2
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{4}}, LongRNA{2}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{2}}, LongDNA{4}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{RNAAlphabet{4}}, LongDNA{4}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer{DNAAlphabet{2}}, LongRNA{4}(random_rna_kmer(len))) for _ in 1:reps])

            # BioSequence Construction
            #   Check that kmers can be constructed from a BioSequence
            #   BioSequence → Kmer → BioSequence
            @test all(Bool[check_biosequence_construction(Kmer{DNAAlphabet{2},len}, LongSequence{DNAAlphabet{2}}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{DNAAlphabet{4},len}, LongSequence{DNAAlphabet{4}}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{DNAAlphabet{2},len}, LongSequence{DNAAlphabet{4}}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{DNAAlphabet{4},len}, LongSequence{DNAAlphabet{2}}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{RNAAlphabet{2},len}, LongSequence{RNAAlphabet{2}}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{RNAAlphabet{4},len}, LongSequence{RNAAlphabet{4}}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{RNAAlphabet{2},len}, LongSequence{RNAAlphabet{4}}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer{RNAAlphabet{4},len}, LongSequence{RNAAlphabet{2}}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(AAKmer{len},              LongSequence{AminoAcidAlphabet}(random_aa(len)))    for _ in 1:reps])

            # Check Kmer(::BioSequence) construction
            @test all(Bool[check_longsequence_construction(Kmer,                     LongRNA{4}(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer,                     LongDNA{2}(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_longsequence_construction(Kmer,                     LongAA(random_rna_kmer(len))) for _ in 1:reps])

            # Construction from element arrays
            #   Check that kmers can be constructed from an array of elements
            #   Vector{T} → Kmer{A,K,N} → Vector{T}
            @test all(Bool[check_nucarray_kmer(Kmer{DNAAlphabet{2},len}, random_dna_symbols(len, [0.25, 0.25, 0.25, 0.25, 0.0])) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Kmer{DNAAlphabet{4},len}, random_dna_symbols(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Kmer{RNAAlphabet{2},len}, random_rna_symbols(len, [0.25, 0.25, 0.25, 0.25, 0.0])) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Kmer{RNAAlphabet{4},len}, random_rna_symbols(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(AAKmer{len},              random_aa_symbols(len))     for _ in 1:reps])
            
            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(Kmer{DNAAlphabet{2},len}, DNAAlphabet{2},    random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{DNAAlphabet{4},len}, DNAAlphabet{4},    random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{DNAAlphabet{2},len}, DNAAlphabet{4},    random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{DNAAlphabet{4},len}, DNAAlphabet{2},    random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{RNAAlphabet{2},len}, RNAAlphabet{2},    random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{RNAAlphabet{4},len}, RNAAlphabet{4},    random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{RNAAlphabet{2},len}, RNAAlphabet{4},    random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer{RNAAlphabet{4},len}, RNAAlphabet{2},    random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(AAKmer{len},              AminoAcidAlphabet, random_aa(len))       for _ in 1:reps])
        end
    end

    @test_throws MethodError Kmer() # can't construct 0-mer using `Kmer()`
    @test_throws ArgumentError DNAKmer(dna"") # 0-mers not allowed
    @test_throws ArgumentError AAKmer(aa"") # 0-mers not allowed
    @test_throws ArgumentError DNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError AAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws BioSequences.EncodeError Kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws BioSequences.EncodeError Kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws BioSequences.EncodeError RNAKmer(rna"ACGNU") # no Ns in 2-bit nucleic acid kmers
    @test_throws BioSequences.EncodeError DNAKmer(dna"ACGNT") # no Ns in 2-bit nucleic acid kmers
    @test_throws MethodError Kmer(RNA_A, DNA_A) # no mixing of RNA and DNA

    @testset "From strings" begin
        @test DNAKmer("ACTG") == DNAKmer(LongDNA{4}("ACTG"))
        @test RNAKmer("ACUG") == RNAKmer(LongRNA{4}("ACUG"))

        # N is not allowed in Kmers
        @test_throws Exception DNAMmer("ACGTNACGT")
        @test_throws Exception RNAKmer("ACGUNACGU")

        # Test string literals
        @test mer"ACTG"dna == DNAKmer(LongDNA{4}("ACTG"))
        @test mer"AVBM"aa  == AAKmer(LongAA("AVBM"))
        @test isa(mer"ACGT"dna, DNAKmer{4})
        @test isa(mer"AVBM"aa,  AAKmer{4})
        @test_throws LoadError eval(:(mer"ACGN"dna))
        @test_throws LoadError eval(:(mer"ACG-"dna))
    end
    
    @testset "Capacity" begin
        @test Kmers.capacity(DNAKmer(random_dna_kmer(10))) == 32
        @test Kmers.capacity(RNAKmer(random_rna_kmer(10))) == 32
        @test Kmers.capacity(DNAKmer(random_dna_kmer(32))) == 32
        @test Kmers.capacity(RNAKmer(random_rna_kmer(32))) == 32
        @test Kmers.capacity(DNAKmer(random_dna_kmer(33))) == 64
        @test Kmers.capacity(AAKmer(random_aa(8))) == 8
        @test Kmers.capacity(AAKmer(random_aa(10))) == 16
    end
    
    @testset "N unused" begin
        @test Kmers.n_unused(DNAKmer(random_dna_kmer(10))) == 22
        @test Kmers.n_unused(RNAKmer(random_rna_kmer(10))) == 22
        @test Kmers.n_unused(DNAKmer(random_dna_kmer(32))) == 0
        @test Kmers.n_unused(RNAKmer(random_rna_kmer(32))) == 0
        @test Kmers.n_unused(DNAKmer(random_dna_kmer(33))) == 31
        @test Kmers.n_unused(AAKmer(random_aa(8))) == 0
        @test Kmers.n_unused(AAKmer(random_aa(10))) == 6
    end
end
