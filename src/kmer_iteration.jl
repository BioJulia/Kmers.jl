

#=
@inline function Base.iterate(it::SpacedKmerIterator{Kmer{A,K,N},S}, state=(it.start - it.increment, 1, 0, blank_ntuple(Kmer{A,K,N}), blank_ntuple(Kmer{A,K,N}))
    ) where {A,K,N,S<:LongSequence{<:NucleicAcidAlphabet{4}}}
    i, pos, filled, fwkmer, rvkmer = state
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(UInt8, @inbounds getindex(it.seq, i))
        @inbounds fbits = UInt64(kmerbits[nt + 1])
        rbits = ~fbits & typeof(fbits)(0x03)
        if fbits == 0xff # ambiguous
            filled = 0
            # Find the beginning of next possible kmer after i
            pos = i + it.step - Core.Intrinsics.urem_int(i - pos, it.step)
            i = pos - 1
        else
            filled += 1
            fwkmer = leftshift_carry(fwkmer, 2, fbits)
            rvkmer = rightshift_carry(rvkmer, 2, UInt64(rbits) << (62 - (64N - 2K)))
        end
        if filled == K
            state = (i, i - K + 1 + it.step, it.filled, fwkmer, rvkmer)
            return KmerAt(pos, Kmer{A,K,N}(fwkmer), Kmer{A,K,N}(rvkmer)), state
        end
        i += 1
    end
    return nothing
end
=#