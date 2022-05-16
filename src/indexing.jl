@inline BioSequences.encoded_data_eltype(::Type{<:Kmer}) = UInt64

@inline function BioSequences.extract_encoded_element(seq::Kmer, i::Integer)
    bi = BioSequences.bitindex(seq, i % UInt)
    return BioSequences.extract_encoded_element(bi, seq.data)
end

@inline Base.copy(seq::Kmer) = typeof(seq)(seq.data)

@inline encoded_data(x::Kmer) = x.data

@inline BioSequences.bitindex(seq::Kmer, i::Integer) = BioSequences.bitindex(BioSequences.BitsPerSymbol(seq), BioSequences.encoded_data_eltype(typeof(seq)), i + n_unused(seq))

