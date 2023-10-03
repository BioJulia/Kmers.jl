"""
    AbstractKmerIterator{A <: Alphabet, K}

Iterates `Kmer{A, K}`.
Functions to implement:
* `Base.iterate`

Optional functions:
* `source_type`
* `load_source`
"""
abstract type AbstractKmerIterator{A <: Alphabet, K} end

function Base.eltype(::Type{<:AbstractKmerIterator{A, K}}) where {A, K}
    Kmer{A, K, n_coding_elements(Kmer{A, K})}
end

"""
    source_type(::Type{<:AbstractKmerIterator})::Type

Get the type of the data source that kmers are extracted from
"""
function source_type end

"""
    load_source(x::AbstractKmerIterator)::source_type(typeof(x))

Get the data source from the kmer iterator.
"""
function load_source end

function used_source(R::RecodingScheme, s)
    if R isa AsciiEncode && s isa Union{String, SubString{String}}
        codeunits(s)
    else
        s
    end
end

@noinline throw_bad_byte_error(b::UInt8) =
    error("Cannot interpret byte $(repr(b)) as nucleotide")

const ASCII_SKIPPING_LUT = let
    v = fill(0xff, 256)
    for (i, s) in [(0, "Aa"), (1, "cC"), (2, "gG"), (3, "TtUu")], c in s
        v[UInt8(c) + 1] = i
    end
    for c in "-MRSVWYHKDBN"
        v[UInt8(c) + 1] = 0xf0
        v[UInt8(lowercase(c)) + 1] = 0xf0
    end
    Tuple(v)
end

const FOURBIT_COMPLEMENT_LUT = let
    v = fill(0x00, 16)
    for i in alphabet(DNA)
        v[reinterpret(UInt8, i) + 0x01] = reinterpret(UInt8, complement(i))
    end
    Tuple(v)
end
