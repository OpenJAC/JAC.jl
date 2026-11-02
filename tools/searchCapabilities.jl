#
# Search CAPABILITIES.md for a term and print the matching entries, each with its "###" heading, its full
# body and its Keywords line, rather than a bare grep line with no context.
#
# Usage:   julia tools/searchCapabilities.jl "<term>"
#          julia tools/searchCapabilities.jl hyperfine
#
if  length(ARGS) != 1
    println("Usage:  julia tools/searchCapabilities.jl \"<term>\"")
    exit(1)
end

term = lowercase(ARGS[1])
path = joinpath(@__DIR__, "..", "CAPABILITIES.md")
text = read(path, String)

# Split the file into entries at each "### " heading; the first chunk (before any heading) is the preamble
# and is never a match target.
chunks = split(text, r"(?=^### )"m)

function  main()
    found = false
    for  chunk in chunks[2:end]
        if  occursin(term, lowercase(chunk))
            heading = split(chunk, '\n')[1]
            println("="^100)
            println(heading)
            println("-"^100)
            print(chunk)
            println()
            found = true
        end
    end

    return( found )
end

if  !main()
    println("No entry in CAPABILITIES.md matches \"$(ARGS[1])\".")
    println("Try a shorter or more general term -- e.g. the physics quantity rather than a JAC type name.")
end
