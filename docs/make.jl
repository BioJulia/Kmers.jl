using Documenter, Kmers

DocMeta.setdocmeta!(
    Kmers,
    :DocTestSetup,
    :(using BioSequences, Kmers, Test);
    recursive=true,
)

makedocs(;
    modules=[Kmers],
    format=Documenter.HTML(),
    sitename="Kmers.jl",
    pages=[
        "Home" => "index.md",
        "The Kmer type" => "kmers.md",
        "Translation" => "translation.md",
        # The kmer type (construction, indexing)
        # Kmer iteration
        # Translation (revtrans also)
        # FAQ (why not compare to bioseq, why no unambig canonical)
    ],
    authors="Jakob Nybo Nissen, Sabrina J. Ward, The BioJulia Organisation and other contributors.",
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/BioJulia/Kmers.jl.git",
    push_preview=true,
    deps=nothing,
    make=nothing,
)
