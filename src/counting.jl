function BioSequences._n_gc(x::Kmer{<:TwoBit})
    mask = 0x5555555555555555 % UInt
    n = 0
    for i in BioSequences.encoded_data(x)
        n += count_ones((i ⊻ (i >>> 1)) & mask)
    end
    return n
end

@inline function BioSequences.count_symbol(x::Oligo{A, U}, sym::BioSymbol) where {A, U}
    iszero(BioSequences.bits_per_symbol(x)) && return length(x)
    iszero(x.x) && return 0
    enc = (BioSequences.encode(Alphabet(x), sym)) % UInt

    # Extend sub-word values to a machine word, so they take the same fast path
    # as a one-word Oligo. The mask must be extended too: when counting the
    # zero encoding, its complement prevents the newly introduced high bits from
    # being counted as padding.
    if sizeof(U) <= sizeof(UInt)
        u = x.x % UInt
        mask = coding_mask(x) % UInt
        u = if iszero(enc)
            u | ~mask
        else
            u & mask
        end
        return BioSequences.count_encoding(u, enc, BioSequences.BitsPerSymbol(x))
    end

    u = x.x
    mask = coding_mask(x)
    # If encoding is zeroed, make sure to set all noncoding bits to 1,
    # so we don't count those.
    u = if iszero(enc)
        u | ~mask
    else
        u & mask
    end
    count = 0
    for i in 0:(div(sizeof(x), sizeof(UInt)) - 1)
        uuint = right_shift(u, i * 8 * sizeof(UInt)) % UInt
        count += BioSequences.count_encoding(uuint, enc, BioSequences.BitsPerSymbol(x))
    end
    return count
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
