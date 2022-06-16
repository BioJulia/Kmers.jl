using Documenter, Kmers

makedocs(
    format = Documenter.HTML(),
    sitename = "Kmers.jl",
    pages = [
        "Home"                           => "index.md",
        "Kmer types"                     => "kmer_types.md",
        "Constructing kmers"             => "construction.md",
        "Indexing & modifying kmers"     => "transforms.md",
        "Predicates"                     => "predicates.md",
        "Random kmers"                   => "random.md",
        "Iterating over Kmers"           => "iteration.md",
        "Translation"                    => "translate.md",
        #"Pattern matching and searching" => "sequence_search.md",
        #"Iteration"                      => "iteration.md",
        #"Counting"                       => "counting.md",
        #"I/O"                            => "io.md",
        #"Interfaces"                     => "interfaces.md"
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/Kmers.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)
