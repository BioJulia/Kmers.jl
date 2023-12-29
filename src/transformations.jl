function Base.reverse(x::Kmer)
    # ( ABC, DEFG) # reverse each element
    # (CBA , GFED) # reverse elements
    # (GFED, CBA ) # rightshift carry a zero
    # ( GFE, DBCA) # final result
    Bps = BioSequences.BitsPerSymbol(Alphabet(x))
    data = map(i -> BioSequences.reversebits(i, Bps), reverse(x.data))
    (_, data) = rightshift_carry(data, bits_unused(typeof(x)), zero(UInt))
    typeof(x)(unsafe, data)
end

# For this method, we don't need to mask the unused bits, because the complement of
# 0x0 == DNA_Gap is still DNA_Gap 
function BioSequences.complement(x::Kmer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}})
    isempty(x) && return x
    data = map(i -> BioSequences.complement_bitpar(i, Alphabet(x)), x.data)
    typeof(x)(unsafe, data)
end

# For this method we do need to mask unused bits, unlike above
function BioSequences.complement(x::Kmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    isempty(x) && return x
    data = map(i -> BioSequences.complement_bitpar(i, Alphabet(x)), x.data)
    (head, tail...) = data
    typeof(x)(unsafe, ((head & get_mask(typeof(x))), tail...))
end

# Generic fallback
function BioSequences.complement(x::Kmer{<:NucleicAcidAlphabet})
    typeof(x)((complement(i) for i in x))
end

# TODO: Should this be the generic BioSequence def in BioSequences.jl?
function BioSequences.reverse_complement(x::Kmer)
    @inline(reverse(@inline(complement(x))))
end

function BioSequences.canonical(x::Kmer)
    rc = reverse_complement(x)
    ifelse(x < rc, x, rc)
end

BioSequences.iscanonical(x::Kmer) = x <= reverse_complement(x)

function BioSequences.translate(
    seq::Kmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}};
    code::BioSequences.GeneticCode=BioSequences.standard_genetic_code,
    allow_ambiguous_codons::Bool=true, # noop in this method
    alternative_start::Bool=false,
)
    n_aa, remainder = divrem(length(seq), 3)
    iszero(remainder) ||
        error("LongRNA length is not divisible by three. Cannot translate.")
    N = n_coding_elements(Kmer{AminoAcidAlphabet, n_aa})
    T = Kmer{AminoAcidAlphabet, n_aa, N}
    data = zero_tuple(T)
    @inbounds for i in 1:n_aa
        a = seq[3i - 2]
        b = seq[3i - 1]
        c = seq[3i - 0]
        codon = BioSequences.unambiguous_codon(a, b, c)
        aa = code[codon]
        carry = UInt(reinterpret(UInt8, aa))
        (_, data) =
            leftshift_carry(data, BioSequences.bits_per_symbol(AminoAcidAlphabet()), carry)
    end
    result = T(unsafe, data)
    if alternative_start && !iszero(ksize(typeof(seq)))
        return setindex(result, 1, AA_M)
    else
        return result
    end
end

function BioSequences.translate(
    seq::Kmer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}};
    code::BioSequences.GeneticCode=BioSequences.standard_genetic_code,
    allow_ambiguous_codons::Bool=true, # noop in this method
    alternative_start::Bool=false,
)
    n_aa, remainder = divrem(length(seq), 3)
    iszero(remainder) ||
        error("LongRNA length is not divisible by three. Cannot translate.")
    N = n_coding_elements(Kmer{AminoAcidAlphabet, n_aa})
    T = Kmer{AminoAcidAlphabet, n_aa, N}
    data = zero_tuple(T)
    @inbounds for i in 1:n_aa
        a = reinterpret(RNA, seq[3i - 2])
        b = reinterpret(RNA, seq[3i - 1])
        c = reinterpret(RNA, seq[3i - 0])
        aa = if isgap(a) | isgap(b) | isgap(c)
            error("Cannot translate nucleotide sequences with gaps.")
        elseif iscertain(a) & iscertain(b) & iscertain(c)
            code[BioSequences.unambiguous_codon(a, b, c)]
        else
            BioSequences.try_translate_ambiguous_codon(code, a, b, c, allow_ambiguous_codons)
        end
        carry = UInt(reinterpret(UInt8, aa))
        (_, data) =
            leftshift_carry(data, BioSequences.bits_per_symbol(AminoAcidAlphabet()), carry)
    end
    result = T(unsafe, data)
    if alternative_start && !iszero(ksize(typeof(seq)))
        return setindex(result, 1, AA_M)
    else
        return result
    end
end
