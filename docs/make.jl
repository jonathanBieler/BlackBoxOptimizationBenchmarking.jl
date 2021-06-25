using BlackBoxOptimizationBenchmarking
using Documenter

DocMeta.setdocmeta!(BlackBoxOptimizationBenchmarking, :DocTestSetup, :(using BlackBoxOptimizationBenchmarking); recursive=true)

makedocs(;
    modules=[BlackBoxOptimizationBenchmarking],
    authors="Jonathan Bieler <jonathan.bieler@alumni.epfl.ch> and contributors",
    repo="https://github.com/jonathanBieler/BlackBoxOptimizationBenchmarking.jl/blob/{commit}{path}#{line}",
    sitename="BlackBoxOptimizationBenchmarking.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jonathanBieler.github.io/BlackBoxOptimizationBenchmarking.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jonathanBieler/BlackBoxOptimizationBenchmarking.jl",
    devbranch = "main",
)
