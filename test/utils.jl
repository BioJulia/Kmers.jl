# Return a random DNA/RNA sequence of the given length.
function random_seq(n::Integer, nts, probs, outtype = String)
    cumprobs = cumsum(probs)
    x = Vector{Char}(undef, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return outtype(x)
end

function random_seq(::Type{A}, n::Integer) where {A<:Alphabet}
    # TODO: Resolve the use of symbols(A()).
    nts = symbols(A())
    probs = Vector{Float64}(undef, length(nts))
    fill!(probs, 1 / length(nts))
    return LongSequence{A}(random_seq(n, nts, probs))
end

function random_dna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs)
end

function random_rna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs)
end

function random_aa(len)
    return random_seq(len,
        ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ],
        push!(fill(0.049, 20), 0.02))
end

function random_dna_symbols(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs, Vector{DNA})
end

function random_rna_symbols(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs, Vector{RNA})
end

function random_rna_symbols(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs, Vector{RNA})
end

function random_aa_symbols(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
     'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ], probs, Vector{AminoAcid})
end

function random_dna_kmer(len)
    return random_dna(len, [0.25, 0.25, 0.25, 0.25])
end

function random_rna_kmer(len)
    return random_rna(len, [0.25, 0.25, 0.25, 0.25])
end

function dna_complement(seq::AbstractString)
    seqc = Vector{Char}(undef, length(seq))
    complementer = Dict(zip("-ACGTSWYRKMDVHBN", "-TGCASWRYMKHBDVN"))
    for (i, c) in enumerate(seq)
        seqc[i] = complementer[c]
    end
    return String(seqc)
end

function rna_complement(seq::AbstractString)
    seqc = Vector{Char}(undef, length(seq))
    complementer = Dict(zip("-ACGUSWYRKMDVHBN", "-UGCASWRYMKHBDVN"))
    for (i, c) in enumerate(seq)
        seqc[i] = complementer[c]
    end
    return String(seqc)
end