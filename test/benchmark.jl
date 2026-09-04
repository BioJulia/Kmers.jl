module TestBenchmarks

using Kmers

extract_kmer(::Type{<:Kmers.AbstractKmerIterator}, x) = x
extract_kmer(::Type{<:FwRvIterator}, x) = first(x)
extract_kmer(::Type{<:UnambiguousKmers}, x) = first(x)

function reducer(it)
    y = 0
    for i in it
        y ⊻= extract_kmer(typeof(it), i).data[1]
    end
    return y
end

using Random: Xoshiro
using BioSequences
rng = Xoshiro(439824)

const N = 10_000_000

const seq_2bit = randseq(DNAAlphabet{2}(), N)
const seq_4bit = randrnaseq(N)
const seq_4bit_unambiguous = LongDNA{4}(seq_2bit)
const seq_aa = randaaseq(N)
const data = String(rand(codeunits("AaCcGgTt"), N))

const names = [
    ("2-bit LongSequence", seq_2bit),
    ("4-bit LongSequence", seq_4bit),
    ("AA LongSequence", seq_aa),
    ("String", data),
]

println("Time to iterate over $N symbols:\n")
println("FwKmers")
for (name, seq) in names[1:3]
    it = FwKmers{typeof(Alphabet(seq)), 7}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = FwDNAMers{7}(data)
y = reducer(it)
@time "    String" reducer(it)

println("FwRvIterator")
for (name, seq) in names[1:2]
    it = FwRvIterator{typeof(Alphabet(seq)), 7}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = FwRvIterator{DNAAlphabet{2}, 7}(data)
y = reducer(it)
@time "    String" reducer(it)

println("CanonicalKmers")
for (name, seq) in names[1:2]
    it = CanonicalKmers{typeof(Alphabet(seq)), 7}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = CanonicalDNAMers{7}(data)
y = reducer(it)
@time "    String" reducer(it)

println("UnambiguousKmers")
for (name, seq) in names[1:2]
    it = UnambiguousRNAMers{7}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = UnambiguousRNAMers{7}(data)
y = reducer(it)
@time "    String" reducer(it)

println("SpacedKmers, step 5")
for (name, seq) in names[1:3]
    it = SpacedKmers{typeof(Alphabet(seq)), 7, 5}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = SpacedDNAMers{7, 5}(data)
y = reducer(it)
@time "    String" reducer(it)

println("SpacedKmers, step 7")
for (name, seq) in names[1:3]
    it = SpacedKmers{typeof(Alphabet(seq)), 7, 7}(seq)
    y = reducer(it)
    @time "    " * name reducer(it)
end
it = SpacedDNAMers{7, 7}(data)
y = reducer(it)
@time "    String" reducer(it)

function unsafe_extract_minimizer(seq::LongDNA{2}, i::Int, ::Val{K}, ::Val{W}) where {K, W}
    T = derive_type(Kmer{DNAAlphabet{2}, K})
    kmer = Kmers.unsafe_extract(Kmers.Copyable(), T, seq, i)
    hash = fx_hash(kmer)
    for offset in 0:(W - 2)
        new_kmer =
            Kmers.unsafe_shift_from(Kmers.Copyable(), kmer, seq, i + K + offset, Val(1))
        new_hash = fx_hash(new_kmer)
        if new_hash < hash
            hash = new_hash
            kmer = new_kmer
        end
    end
    return kmer
end

function benchmark_minimizer(seq)
    y = 0
    for i in 1:20:(length(seq) - 19)
        mer = unsafe_extract_minimizer(seq, i, Val{8}(), Val{20}())
        y ⊻= mer.data[1]
    end
    return y
end

println("Minimizer")
benchmark_minimizer(seq_2bit)
@time "    2-bit LongSequence" benchmark_minimizer(seq_2bit)

function benchmark_extract(recoding, ::Type{T}, seq, from) where {T}
    h = zero(UInt)
    stride = BioSequences.symbols_per_data_element(seq)
    for i in 0:249_999
        h ⊻= fx_hash(Kmers.unsafe_extract(recoding, T, seq, from + i * stride))
    end
    return h
end

println("unsafe_extract")
for (name, recoding, T, seq, from) in [
        ("Copyable, aligned small", Kmers.Copyable(), DNAKmer{16, 1}, seq_2bit, 1),
        ("Copyable, misaligned small", Kmers.Copyable(), DNAKmer{16, 1}, seq_2bit, 2),
        ("Copyable, aligned multiword", Kmers.Copyable(), DNAKmer{64, 2}, seq_2bit, 1),
        ("Copyable, misaligned multiword", Kmers.Copyable(), DNAKmer{64, 2}, seq_2bit, 2),
        (
            "TwoToFour, aligned small",
            Kmers.TwoToFour(),
            Kmer{DNAAlphabet{4}, 16, 1},
            seq_2bit,
            1,
        ),
        (
            "TwoToFour, misaligned multiword",
            Kmers.TwoToFour(),
            Kmer{DNAAlphabet{4}, 32, 2},
            seq_2bit,
            2,
        ),
        (
            "FourToTwo, aligned small",
            Kmers.FourToTwo(),
            DNAKmer{16, 1},
            seq_4bit_unambiguous,
            1,
        ),
        (
            "FourToTwo, misaligned multiword",
            Kmers.FourToTwo(),
            DNAKmer{64, 2},
            seq_4bit_unambiguous,
            2,
        ),
    ]
    benchmark_extract(recoding, T, seq, from)
    @time "    " * name benchmark_extract(recoding, T, seq, from)
end

end # module
