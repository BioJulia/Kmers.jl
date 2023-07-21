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

# For this method we do
function BioSequences.complement(x::Kmer{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    isempty(x) && return x
    data = map(i -> BioSequences.complement_bitpar(i, Alphabet(x)), x.data)
    (head, tail...) = data
    typeof(x)(unsafe, ((head & get_mask(typeof(x))), tail...))
end

# Generic fallback
function BioSequences.complement(x::Kmer{<:NucleicAcidAlphabet})
    construct_generic_unchecked(Base.HasLength(), typeof(x), (complement(i) for i in x))
end

# TODO: Should this be the generic BioSequence def in BioSequences.jl?
function BioSequences.reverse_complement(x::Kmer)
    reverse(complement(x))
end

function BioSequences.canonical(x::Kmer)
    rc = reverse_complement(x)
    ifelse(x < rc, x, rc)
end

BioSequences.iscanonical(x::Kmer) = x <= reverse_complement(x)

# TODO: Translation