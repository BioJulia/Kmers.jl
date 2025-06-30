function BioSequences._n_gc(x::Kmer{<:TwoBit})
    mask = 0x5555555555555555 % UInt
    n = 0
    for i in BioSequences.encoded_data(x)
        n += count_ones((i ⊻ (i >>> 1)) & mask)
    end
    return n
end

@inline function BioSequences.count_symbol(x::DynamicKmer, sym::BioSymbol)
    iszero(BioSequences.bits_per_symbol(x)) && return length(x)
    u = x.x
    iszero(u) && return 0
    U = utype(typeof(x))
    enc = U(BioSequences.encode(Alphabet(x), sym))
    mask = coding_mask(x)
    u = if iszero(enc)
        u | ~mask
    else
        u & mask
    end
    return BioSequences.count_encoding(u, enc, BioSequences.BitsPerSymbol(x))
end

@inline function BioSequences.count_symbol(x::Kmer, sym::BioSymbol)
    iszero(BioSequences.bits_per_symbol(x)) && return length(x)
    isempty(x) && return 0
    enc = UInt64(BioSequences.encode(Alphabet(x), sym))
    (head, rest...) = x.data
    mask = get_mask(typeof(x))
    head = if iszero(enc)
        head | ~mask
    else
        head
    end
    BPS = BioSequences.BitsPerSymbol(x)
    result = BioSequences.count_encoding(head, enc, BPS)
    for i in rest
        result += BioSequences.count_encoding(i, enc, BPS)
    end
    return result
end
