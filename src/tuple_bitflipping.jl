# These compile to raw CPU instructions and are therefore more
# efficient than simply using << and >>>
@inline function left_shift(x::Unsigned, n::Integer)
    return x << (n & ((sizeof(x) * 8) - 1))
end

@inline function right_shift(x::Unsigned, n::Integer)
    return x >>> (n & ((sizeof(x) * 8) - 1))
end

# When the UInt is shifted n bits, these are the bits
# that are shifted away (carried over)
@inline function left_carry(x::Unsigned, n::Integer)
    return right_shift(x, 8 * sizeof(x) - n)
end

@inline function right_carry(x::Unsigned, n::Integer)
    return left_shift(x, 8 * sizeof(x) - n)
end

# Shift a tuple left nbits, carry over bits between tuple elements, and OR
# the `carry` argument to the right side of the resulting tuple.
# Returns (new_carry, new_tuple)
@inline function leftshift_carry(
        x::Tuple{Vararg{T}},
        nbits::Integer,
        carry::T,
    ) where {T <: Unsigned}
    isempty(x) && return x
    (new_carry, new_tail) = leftshift_carry(tail(x), nbits, carry)
    new_head = left_shift(first(x), nbits) | new_carry
    return (left_carry(first(x), nbits), (new_head, new_tail...))
end

@inline function rightshift_carry(
        x::Tuple{Vararg{T}},
        nbits::Integer,
        carry::T,
    ) where {T <: Unsigned}
    isempty(x) && return x
    new_head = right_shift(first(x), nbits) | right_carry(carry, nbits)
    mask = left_shift(UInt(1), nbits) - 1
    tail_carry = first(x) & mask
    (new_carry, new_tail) = rightshift_carry(tail(x), nbits, tail_carry)
    return (new_carry, (new_head, new_tail...))
end

# Recusion terminator for above
@inline leftshift_carry(::Tuple{}, nbits::Integer, carry::Unsigned) = (carry, ())
@inline rightshift_carry(::Tuple{}, nbits::Integer, carry::Unsigned) = (carry, ())
