using MarkovOptim
using Documenter

DocMeta.setdocmeta!(MarkovOptim, :DocTestSetup, :(using MarkovOptim); recursive=true)

makedocs(;
    modules=[MarkovOptim],
    authors="Violet Hart",
    repo="https://github.com/pjh6654/MarkovOptim.jl/blob/{commit}{path}#{line}",
    sitename="MarkovOptim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pjh6654.github.io/MarkovOptim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pjh6654/MarkovOptim.jl",
)
