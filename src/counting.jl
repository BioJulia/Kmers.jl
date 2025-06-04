function BioSequences._n_gc(x::Kmer{<:TwoBit})
    mask = 0x5555555555555555 % UInt
    n = 0
    for i in BioSequences.encoded_data(x)
        n += count_ones((i ⊻ (i >>> 1)) & mask)
    end
    return n
end
