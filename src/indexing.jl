@inline function BioSequences.extract_encoded_element(seq::Kmer{A}, i::Integer) where A
    T = typeof(seq)
    bps = BioSequences.bits_per_symbol(A()) % UInt
    index = div((i + n_unused(T) - 1) % UInt, per_word_capacity(T) % UInt) + 1
    offset = mod(((elements_in_head(T) - i) * bps) % UInt, 8 * sizeof(UInt))
    mask = UInt(1) << bps - 1
    right_shift(@inbounds(seq.data[index]), offset) & mask
end

# TODO: Index with range, index with bitvector