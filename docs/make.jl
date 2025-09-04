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
        "Exported Functions" => "index.md",
        "Internals" => "internals.md",
    ],
)

deploydocs(;
    repo="github.com/phlaster/SeqFold.jl",
    devbranch="master",
)
