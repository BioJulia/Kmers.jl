module BitIntegersExt

using BitIntegers
using Kmers: Kmers

@inline Base.@constprop :aggressive Base.@assume_effects :foldable function Kmers.get_large_bitsize(
        ::Val{N}
    ) where {N}
    if N < 32
        UInt256
    elseif N < 64
        UInt512
    elseif N < 128
        UInt1024
    else
        error("Kmers does not support BitIntegers larger than UInt1024")
    end
end

Kmers.widen_bitint(::Type{UInt128}) = UInt256
Kmers.widen_bitint(::Type{UInt256}) = UInt512
Kmers.widen_bitint(::Type{UInt512}) = UInt1024

end # module
