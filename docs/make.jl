# push!(LOAD_PATH,"../src/")

using Documenter, JAC

makedocs(;
    modules=[JAC],
    format=Documenter.HTML(prettyurls=false, 
                           size_threshold = 300 * 1024,),
    pages=[
        "Home"      => "index.md",
        "Examples"  => "examples.md",
        "API Reference"   => "api.md",
        "Publications Reference" => "reference.md",
        "News"        => "news.md",
        "Contributors" => "contributors.md",
        "License"   => "license.md",
    ],
    repo="https://github.com/AlokaSahoo/JAC.jl/",
    sitename="JAC.jl",
    authors="Stephan Fritzsche",
    warnonly = Documenter.except(),
)

deploydocs(;
    repo = "github.com/AlokaSahoo/JAC.jl"
)