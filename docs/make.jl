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
        sidebar_sitename=false,
    ),
    checkdocs = :public,
    pages=[
        "Home" => "index.md",
        "Melting Temperature" => "tm.md",
        "Sequence Folding" => "fold.md",
        "Other Functions" => "utils.md",
        "For Python users" => "juliacall.md",
        "Citations" => "citations.md"
    ],
)

deploydocs(;
    repo="github.com/phlaster/SeqFold.jl",
    devbranch="master",
)