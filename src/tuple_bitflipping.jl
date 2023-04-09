
# TODO: this should end up in BioSequences.jl?

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function BioSequences.extract_encoded_element(bidx::BioSequences.BitIndex{N,W}, data::NTuple{n,W}) where {N,n,W}
    @inbounds chunk = data[BioSequences.index(bidx)]
    offchunk = chunk >> (BioSequences.bitwidth(W) - N - BioSequences.offset(bidx))
    return offchunk & BioSequences.bitmask(bidx)
end



"""
    _cliphead(by::Integer, head::UInt64, tail...)

A method used to mask the first `by` MSB's in `head`, before catting it with
tail to return a NTuple.

This is used internally to mask the first `by` bits in the first word of a 
NTuple of UInt64's.

Notably it's used when constructing a Kmer from an existing NTuple of UInt64
"""
@inline function _cliphead(by::Integer, head::UInt64, tail...)
    return (head & (typemax(UInt64) >> by), tail...)
end


"""
    _drophead(x::NTuple{N,UInt64}, drop::Val{D}) where {N,D}

A function used to drop the first N elements from an NTuple.

Intended for internal bit-flipping use.

Notably used in cases where expanding a 2-bit encoding to
a 4-bit encoding results in an unused element at the head
of the tuple. 
"""
@inline function _drophead(x::NTuple{N,UInt64}, drop::Val{D}) where {N,D}
    ntuple(Val{N - D}()) do i
        Base.@_inline_meta
        return @inbounds x[D + i]
    end
end

#=
rightshift_carry & leftshift_carry

These methods are micro-optimised (or should be!!!) for shifting the bits in 
an NTuple of unsigned integers, carrying the bits "shifted off" one word 
over to the next word. The carry can also be "seeded" so as other methods like
pushfirst and pushlast can be efficiently implemented without duplication of code
or less efficient implementations that first shift and then insert an element.
=#

@inline function rightshift_carry(x::NTuple{N,UInt64}, nbits::Integer, prevcarry = zero(UInt64)) where {N}
    return _rightshift_carry(nbits, prevcarry, x...)
end

@inline function _rightshift_carry(nbits::Integer, carry::UInt64, head::UInt64, tail...)
    return ((head >> nbits) | carry, _rightshift_carry(nbits, (head & ((one(UInt64) << nbits) - 1)) << (64 - nbits), tail...)...)
end

@inline _rightshift_carry(nbits::Integer, carry::UInt64) = ()

@inline function leftshift_carry(x::NTuple{N,UInt64}, nbits::Integer, prevcarry::UInt64 = zero(UInt64)) where {N}
    _, newbits = _leftshift_carry(nbits, prevcarry, x...)
    return newbits
end

@inline function _leftshift_carry(nbits::Integer, prevcarry::UInt64, head::UInt64, tail...)
    carry, newtail = _leftshift_carry(nbits, prevcarry, tail...)
    return head >> (64 - nbits), ((head << nbits) | carry, newtail...)
end

@inline _leftshift_carry(nbits::Integer, prevcarry::UInt64) = prevcarry, ()

@inline function _reverse(bpe::BioSequences.BitsPerSymbol{N}, head::UInt64, tail...) where {N}
    return (_reverse(bpe, tail...)..., BioSequences.reversebits(head, bpe))
end

@inline _reverse(::BioSequences.BitsPerSymbol{N}) where {N} = ()

#=
@inline function _reverse(f::F, bpe::BioSequences.BitsPerSymbol{N}, head::UInt64, tail...) where {N,F<:Function}
    return (_reverse(f, bpe, tail...)..., f(reversebits(head, bpe)))
end

@inline _reverse(f::F, ::BioSequences.BitsPerSymbol{N}) where {N,F<:Function} = ()
=#

const transcode_2_to_4_nucs = (function ()

    @info "Computing 2 -> 4 bit translation table..."

    fourbitnucs = (UInt32(1), UInt32(2), UInt32(4), UInt32(8))

    fn(bits, lshift) = fourbitnucs[((bits >> lshift) & 0x03) + 1] << (2 * lshift) 

    function translate2to4bits(bits::UInt16)
        fn(bits, 0)  |
        fn(bits, 2)  |
        fn(bits, 4)  |
        fn(bits, 6)  |
        fn(bits, 8)  |
        fn(bits, 10) |
        fn(bits, 12) |
        fn(bits, 14) 
    end

    ntuple((x) -> translate2to4bits(UInt16(x - 1)), 65536)
end)()

@inline function transcode_kernel(bits::UInt64, from::A, to::B) where
    {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4}}

    @inbounds begin
        b = UInt64(transcode_2_to_4_nucs[((bits >>  0) & 0xFFFF) + 1]) <<  0 |
            UInt64(transcode_2_to_4_nucs[((bits >> 16) & 0xFFFF) + 1]) << 32

        a = UInt64(transcode_2_to_4_nucs[((bits >> 32) & 0xFFFF) + 1]) <<  0 |
            UInt64(transcode_2_to_4_nucs[((bits >> 48) & 0xFFFF) + 1]) << 32
    end
    return (a, b)
end

@inline function transcode_bits(x::NTuple{N,UInt64}, from::A, to::B) where
    {N,A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4}}

    return _transcode_bits(from, to, x...)
end

@inline function _transcode_bits(from::A, to::B, head::UInt64, tail...) where
    {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4}}

    return (transcode_kernel(head, from, to)..., _transcode_bits(from, to, tail...)...)
end

@inline function _transcode_bits(from::A, to::B, head::UInt64) where
    {A<:NucleicAcidAlphabet{2},B<:NucleicAcidAlphabet{4}}

    return transcode_kernel(head, from, to)
end
