
# TODO: this should end up in BioSequences.jl?

#=
"Extract the element stored in a packed bitarray referred to by bidx."
@inline function BioSequences.extract_encoded_element(
    bidx::BioSequences.BitIndex{N, W},
    data::NTuple{n, W},
) where {N, n, W}
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

# rightshift_carry & leftshift_carry

# These methods are micro-optimised (or should be!!!) for shifting the bits in 
# an NTuple of unsigned integers, carrying the bits "shifted off" one word 
# over to the next word. The carry can also be "seeded" so as other methods like
# pushfirst and pushlast can be efficiently implemented without duplication of code
# or less efficient implementations that first shift and then insert an element.

@inline function rightshift_carry(
    x::NTuple{N, UInt64},
    nbits::Integer,
    prevcarry=zero(UInt64),
) where {N}
    return _rightshift_carry(nbits, prevcarry, x...)
end

@inline function _rightshift_carry(nbits::Integer, carry::UInt64, head::UInt64, tail...)
    return (
        (head >> nbits) | carry,
        _rightshift_carry(
            nbits,
            (head & ((one(UInt64) << nbits) - 1)) << (64 - nbits),
            tail...,
        )...,
    )
end

@inline _rightshift_carry(nbits::Integer, carry::UInt64) = ()

@inline function leftshift_carry(
    x::NTuple{N, UInt64},
    nbits::Integer,
    prevcarry::UInt64=zero(UInt64),
) where {N}
    _, newbits = _leftshift_carry(nbits, prevcarry, x...)
    return newbits
end

@inline function _leftshift_carry(nbits::Integer, prevcarry::UInt64, head::UInt64, tail...)
    carry, newtail = _leftshift_carry(nbits, prevcarry, tail...)
    return head >> (64 - nbits), ((head << nbits) | carry, newtail...)
end

@inline _leftshift_carry(nbits::Integer, prevcarry::UInt64) = prevcarry, ()

@inline function _reverse(
    bpe::BioSequences.BitsPerSymbol{N},
    head::UInt64,
    tail...,
) where {N}
    return (_reverse(bpe, tail...)..., BioSequences.reversebits(head, bpe))
end

@inline _reverse(::BioSequences.BitsPerSymbol{N}) where {N} = ()
=#

# These compile to raw CPU instructions and are therefore more
# efficient than simply using << and >>>
@inline function left_shift(x::Unsigned, n::Integer)
    x << (n & ((sizeof(x) * 8) - 1))
end

@inline function right_shift(x::Unsigned, n::Integer)
    x >>> (n & ((sizeof(x) * 8) - 1))
end

# When the UInt is shifted n bits, these are the bits
# that are shifted away (carried over)
@inline function left_carry(x::Unsigned, n::Integer)
    right_shift(x, 8 * sizeof(x) - n)
end

@inline function right_carry(x::Unsigned, n::Integer)
    left_shift(x, 8 * sizeof(x) - n)
end

# Shift a tuple left nbits, carry over bits between tuple elements, and OR
# the `carry` argument to the right side of the resulting tuple.
# Returns (new_carry, new_tuple)
@inline function leftshift_carry(x::Tuple{Vararg{T}}, nbits::Integer, carry::T) where {T <: Unsigned}
    head, tail... = x
    (new_carry, new_tail) = leftshift_carry(tail, nbits, carry)
    new_head = left_shift(head, nbits) | new_carry
    (left_carry(head, nbits), (new_head, new_tail...))
end

@inline function rightshift_carry(x::Tuple{Vararg{T}}, nbits::Integer, carry::T) where {T <: Unsigned}
    head, tail... = x
    new_head = right_shift(head, nbits) | right_carry(carry, nbits)
    mask = left_shift(UInt(1), nbits) - 1
    tail_carry = head & mask
    (new_carry, new_tail) = rightshift_carry(tail, nbits, tail_carry)
    (new_carry, (new_head, new_tail...))
end

# Recusion terminator for above
@inline leftshift_carry(::Tuple{}, nbits::Integer, carry::Unsigned) = (carry, ())
@inline rightshift_carry(::Tuple{}, nbits::Integer, carry::Unsigned) = (carry, ())

