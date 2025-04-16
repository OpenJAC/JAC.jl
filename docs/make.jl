# push!(LOAD_PATH,"../src/")

using Documenter, JAC

makedocs(;
    modules=[JAC],
    format=Documenter.HTML(prettyurls=false, 
                           size_threshold = 300 * 1024,),
    pages=[
        "Home"      => "index.md",
        "Examples"  => "examples.md",
        "JAC API Reference"   => "api.md",
        # "Types"     => "types.md",
        # "Functions" => "functions.md",
        "Reference" => "reference.md",
        "Contributors" => "contributors.md",
        "License"   => "license.md",
    ],
    repo="https://github.com/OpenJAC/JAC.jl/",
    sitename="JAC.jl",
    authors="Stephan Fritzsche",
    warnonly = Documenter.except(),
)

deploydocs(;
    repo = "github.com/AlokaSahoo/JAC.jl"
)