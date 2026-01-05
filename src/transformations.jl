function Base.reverse(x::Kmer)
    # ( ABC, DEFG) # reverse each element
    # (CBA , GFED) # reverse elements
    # (GFED, CBA ) # rightshift carry a zero
    # ( GFE, DBCA) # final result
    Bps = BioSequences.BitsPerSymbol(Alphabet(x))
    data = map(i -> BioSequences.reversebits(i, Bps), reverse(x.data))
    (_, data) = rightshift_carry(data, bits_unused(typeof(x)), zero(UInt))
    return typeof(x)(unsafe, data)
end

# For this method, we don't need to mask the unused bits, because the complement of
# 0x0 == DNA_Gap is still DNA_Gap
function BioSequences.complement(x::Kmer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}})
    isempty(x) && return x
    data = map(i -> BioSequences.complement_bitpar(i, Alphabet(x)), x.data)
    return typeof(x)(unsafe, data)
end

# For this method we do need to mask unused bits, unlike above
function BioSequences.complement(x::Kmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    isempty(x) && return x
    data = map(i -> BioSequences.complement_bitpar(i, Alphabet(x)), x.data)
    return typeof(x)(unsafe, ((first(data) & get_mask(typeof(x))), Base.tail(data)...))
end

# Generic fallback
function BioSequences.complement(x::Kmer{<:NucleicAcidAlphabet})
    return typeof(x)((complement(i) for i in x))
end

function BioSequences.reverse_complement(x::Kmer)
    return @inline(reverse(@inline(complement(x))))
end

function BioSequences.canonical(x::Kmer)
    rc = reverse_complement(x)
    return ifelse(x < rc, x, rc)
end

BioSequences.iscanonical(x::Kmer) = x <= reverse_complement(x)

function BioSequences.translate(
        seq::Kmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true, # noop in this method
        alternative_start::Bool = false,
    )
    iszero(ksize(typeof(seq))) && return mer""a
    n_aa, remainder = divrem(length(seq), 3)
    iszero(remainder) ||
        error("LongRNA length is not divisible by three. Cannot translate.")
    N = n_coding_elements(Kmer{AminoAcidAlphabet, n_aa})
    T = Kmer{AminoAcidAlphabet, n_aa, N}
    data = zero_tuple(T)
    # In the next two lines: If alternative_start, we shift in the encoding of M
    # to first place, then we skip the first 3 nucleotides
    (_, data) = leftshift_carry(data, 8, UInt(0x0c) * alternative_start)
    @inbounds for i in (1 + (3 * alternative_start)):n_aa
        a = seq[3i - 2]
        b = seq[3i - 1]
        c = seq[3i - 0]
        codon = BioSequences.unambiguous_codon(a, b, c)
        aa = code[codon]
        carry = UInt(reinterpret(UInt8, aa))
        (_, data) =
            leftshift_carry(data, BioSequences.bits_per_symbol(AminoAcidAlphabet()), carry)
    end
    return T(unsafe, data)
end

function BioSequences.translate(
        seq::Kmer{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}};
        code::BioSequences.GeneticCode = BioSequences.standard_genetic_code,
        allow_ambiguous_codons::Bool = true,
        alternative_start::Bool = false,
    )
    n_aa, remainder = divrem(length(seq), 3)
    iszero(remainder) ||
        error("LongRNA length is not divisible by three. Cannot translate.")
    N = n_coding_elements(Kmer{AminoAcidAlphabet, n_aa})
    T = Kmer{AminoAcidAlphabet, n_aa, N}
    data = zero_tuple(T)
    # In the next two lines: If alternative_start, we shift in the encoding of M
    # to first place, then we skip the first 3 nucleotides
    (_, data) = leftshift_carry(data, 8, UInt(0x0c) * alternative_start)
    @inbounds for i in (1 + (3 * alternative_start)):n_aa
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
    return T(unsafe, data)
end
