using AtomsBuilder
using Documenter

DocMeta.setdocmeta!(AtomsBuilder, :DocTestSetup, :(using AtomsBuilder); recursive=true)

makedocs(;
    modules=[AtomsBuilder],
    authors="Christoph Ortner <christohortner@gmail.com> and contributors",
    sitename="AtomsBuilder.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaMolSim.github.io/AtomsBuilder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaMolSim/AtomsBuilder.jl",
    devbranch="main",
)
