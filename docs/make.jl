using Documenter, Kmers

DocMeta.setdocmeta!(
    Kmers,
    :DocTestSetup,
    :(using BioSequences, Kmers, Test);
    recursive=true,
)

makedocs(;
    modules=[Kmers],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    sitename="Kmers.jl",
    pages=[
        "Home" => "index.md",
        "The Kmer type" => "kmers.md",
        "Iteration" => "iteration.md",
        "Translation" => "translation.md",
        "Hashing" => "hashing.md",
        "Random kmers" => "random.md",
        "K-mer replacements" => "replacements.md",
        "FAQ" => "faq.md",
        "Cookbook" => ["MinHash" => "minhash.md", "Kmer composition" => "composition.md"],
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
