using SeqFold
using Documenter

DocMeta.setdocmeta!(SeqFold, :DocTestSetup, :(using SeqFold); recursive=true)

makedocs(;
    modules=[SeqFold],
    authors="phlaster <phlaster@users.noreply.github.com>",
    sitename="SeqFold.jl",
    format=Documenter.HTML(;
        canonical="https://phlaster.github.io/SeqFold.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/phlaster/SeqFold.jl",
    devbranch="master",
)
