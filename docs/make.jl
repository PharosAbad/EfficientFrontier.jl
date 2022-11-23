using EfficientFrontier
using Documenter

DocMeta.setdocmeta!(EfficientFrontier, :DocTestSetup, :(using EfficientFrontier); recursive=true)

makedocs(;
    modules=[EfficientFrontier],
    authors="Pharos Abad",
    repo="https://github.com/PharosAbad/EfficientFrontier.jl/blob/{commit}{path}#{line}",
    sitename="EfficientFrontier.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://PharosAbad.github.io/EfficientFrontier.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PharosAbad/EfficientFrontier.jl",
    devbranch="main",
)
