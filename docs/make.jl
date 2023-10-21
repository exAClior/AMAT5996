using AMAT5996
using Documenter

DocMeta.setdocmeta!(AMAT5996, :DocTestSetup, :(using AMAT5996); recursive=true)

makedocs(;
    modules=[AMAT5996],
    authors="Yusheng Zhao <yushengzhao2020@outlook.com> and contributors",
    repo="https://github.com/exAClior/AMAT5996.jl/blob/{commit}{path}#{line}",
    sitename="AMAT5996.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://exAClior.github.io/AMAT5996.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/exAClior/AMAT5996.jl",
    devbranch="main",
)
