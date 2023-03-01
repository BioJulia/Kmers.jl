const N_AA = length(AminoAcidAlphabet())

"""
    CodonSet <: AbstractSet{RNACodon}

A small, immutable set of `RNACodon`.

Create an empty set using `CodonSet()`, or from an iterable of `RNACodon`
using `CodonSet(itr)`.
Because `CodonSet` is immutable, use `push` instead of `push!`, and use
the non-mutating set operations `union`, `setdiff`, etc.

# Examples
```jldoctest
julia> v = (mer"UAG"r, mer"GGA"r, mer"UUU"r);

julia> Set(CodonSet(v)) == Set(v)
true

julia> union(CodonSet(v), CodonSet([mer"GAG"r]))
Kmers.CodonSet with 4 elements:
  GAG
  GGA
  UAG
  UUU
```

See also: `push`
"""
struct CodonSet <: AbstractSet{RNACodon}
    x::UInt64

    CodonSet(x::UInt64, ::Unsafe) = new(x)
end
CodonSet() = CodonSet(UInt64(0), Unsafe())
CodonSet(itr) = foldl(push, itr, init=CodonSet())

function Base.iterate(x::CodonSet, s::UInt64=x.x)
    codon = RNACodon((trailing_zeros(s) % UInt64,))
    iszero(s) ? nothing : (codon, s & (s-1))
end

function push(s::CodonSet, x::RNACodon)
    CodonSet(s.x | (UInt64(1) << (x.data[1] & 63)), Unsafe())
end

Base.length(x::CodonSet) = count_ones(x.x)
Base.in(c::RNACodon, s::CodonSet) = isodd(s.x >>> (c.data[1] & 63))
delete(s::CodonSet, x::RNACodon) = setdiff(s, CodonSet((x,)))
Base.issubset(a::CodonSet, b::CodonSet) = isempty(setdiff(a, b))
Base.filter(f, s::CodonSet) = CodonSet(Iterators.filter(f, s))
Base.setdiff(a::CodonSet, b::Vararg{CodonSet}) = CodonSet(a.x & ~(union(b...).x), Unsafe())

for (name, f) in [(:union, |), (:intersect, &), (:symdiff, ⊻)]
    @eval function Base.$(name)(a::CodonSet, b::Vararg{CodonSet}) 
        CodonSet(mapreduce(i -> i.x, $f, b, init=a.x), Unsafe())
    end
end

"""
    ReverseGeneticCode <: AbstractDict{AminoAcid, CodonSet}

A mapping from an amino acid `aa` to the `CodonSet` of all codons that
translate to `aa`. Conceptually, the inverse of a `BioSequences.GeneticCode`.
Used by `reverse_translate`.

`AA_Gap` cannot be translated. Ambiguous amino acids translate to the union
of what their constituent amino acids translate to.
Pyrrolysine and selenocysteine translate to `CodonSet` containing `UAG` and `UGA`,
respectively, whereas they are not translated to in most forward genetic codes.
For these reasons, a the mapping through `ReverseGeneticCode` is not exactly
inverse of the mapping through `GeneticCode`

# Examples
```jldoctest
julia> code = ReverseGeneticCode(BioSequences.candidate_division_sr1_genetic_code);

julia> code[AA_E]
Kmers.CodonSet with 2 elements:
  GAA
  GAG

julia> code[AA_Gap]
ERROR: Cannot reverse translate element: -
[...]
```

See also: [`reverse_translate`](@ref)
"""
struct ReverseGeneticCode <: AbstractDict{AminoAcid, CodonSet}
    name::String
    sets::NTuple{N_AA-1, CodonSet}
end

function ReverseGeneticCode(x::BioSequences.GeneticCode)
    ind(aa::AminoAcid) = reinterpret(UInt8, aa) + 1

    sets = fill(CodonSet(), N_AA-1)
    x_set = CodonSet()
    for i in Int64(0):Int64(63)
        aa = x.tbl[i + 1]
        codon = RNACodon((i % UInt64,))
        sets[ind(aa)] = push(sets[ind(aa)], codon)
        if aa !== AA_Term
            x_set = push(x_set, codon)
        end
    end

    # Ambiguous amino acids
    for (n, (a, b)) in [(AA_B, (AA_D, AA_N)), (AA_J, (AA_I, AA_L)), (AA_Z, (AA_E, AA_Q))]
        sets[ind(n)] = sets[ind(a)] ∪ sets[ind(b)]
    end

    # AA_X codes for all amino acids, except AA_Term
    sets[ind(AA_X)] = x_set

    # Pyrrolysine and selenocysteine are never part of the "forward" genetic
    # code, but can be unambiguously resolved in the reverse genetic code.
    sets[ind(AA_U)] = CodonSet((mer"UGA"r,))
    sets[ind(AA_O)] = CodonSet((mer"UAG"r,))

    ReverseGeneticCode(x.name, Tuple(sets))
end

const rev_standard_genetic_code = ReverseGeneticCode(
    BioSequences.standard_genetic_code
)

function Base.getindex(s::ReverseGeneticCode, a::AminoAcid)
    if reinterpret(UInt8, a) > (N_AA - 2) # cannot translate gap
        error("Cannot reverse translate element: ", a)
    end
    @inbounds s.sets[reinterpret(UInt8, a) + 1]
end

Base.length(c::ReverseGeneticCode) = length(c.sets)
function Base.iterate(c::ReverseGeneticCode, s=1)
    s > length(c.sets) && return nothing
    return (reinterpret(AminoAcid, (s-1)%UInt8) => c.sets[s], s+1)
end

"""
    reverse_translate!(v::Vector{CodonSet}, s::AASeq code=rev_standard_genetic_code)

Reverse-translates `s` under the reverse genetic code `code`, putting the result in `v`.

See also: [`reverse_translate`](@ref)
"""
function reverse_translate!(
    v::Vector{CodonSet},
    seq::AASeq,
    code=rev_standard_genetic_code
)
    resize!(v, length(seq))
    @inbounds for i in eachindex(v)
        v[i] = code[seq[i]]
    end
    v
end

"""
    reverse_translate(s::Union{AminoAcid, AASeq}, code=rev_standard_genetic_code)

Reverse-translates sequence or amino acid `s` under `code::ReverseGeneticCode`
If `s` is an `AminoAcid`, return a `CodonSet`.
If `s` is an `AASeq`, return `Vector{CodonSet}`.

# Examples
```jldoctest
julia> reverse_translate(AA_W)
Kmers.CodonSet with 1 element:
  UGG

julia> v = reverse_translate(aa"MMLVQ");

julia> typeof(v)
Vector{Kmers.CodonSet}

julia> v[4]
Kmers.CodonSet with 4 elements:
  GUA
  GUC
  GUG
  GUU
[...]
```

See also: [`reverse_translate!`](@ref), [`ReverseGeneticCode`](@ref)
"""
function reverse_translate end

function reverse_translate(seq::AASeq, code=rev_standard_genetic_code)
    reverse_translate!(Vector{CodonSet}(undef, length(seq)), seq, code)
end

reverse_translate(aa::AminoAcid, code=rev_standard_genetic_code) = code[aa]
