# push!(LOAD_PATH,"../src/")

using Documenter, JAC

makedocs(;
    modules=[JAC],
    format=Documenter.HTML(prettyurls=false),
    pages=[
        "Home"      => "index.md",
        "Examples"  => "examples.md",
        "Types"     => "types.md",
        "Functions" => "functions.md",
        "Reference" => "reference.md",
        "Contributers" => "contributors.md",
        "License"   => "license.md",
    ],
    repo="https://github.com/OpenJAC/JAC.jl/",
    sitename="JAC.jl",
    authors="Stephan Fritzsche",
    warnonly = Documenter.except(),
)

