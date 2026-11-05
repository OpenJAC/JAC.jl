#
# Re-derive CAPABILITIES.md mechanically from two things that already exist and are already authoritative:
#   1. examples/examples.jl's own one-line description of every example branch;
#   2. the `# Last successful:  DD-Mon-YYYY` dates already recorded in examples/*.jl -- the one date marker
#      in this codebase a HUMAN, not a tool, is allowed to write (CLAUDE.md, Rule 7).
#
# Nothing here invents a date, a module name, or a caveat. Where a fact cannot be derived, the row says so
# plainly (an omitted field, or "unknown/undated") rather than guessing -- a register that oversells an
# unverified branch would be worse than no register at all.
#
# The only hand-maintained input is tools/capabilityNotes.jl -- titles for the few branch-letter groups
# examples.jl itself never names, module-name overrides verified against src/, and a handful of restriction
# hints and keyword synonyms. Everything else below is re-derived every time this script runs.
#
# Usage:   julia --project=. tools/buildCapabilityRegister.jl           (prints to stdout, no file written)
#          julia --project=. tools/buildCapabilityRegister.jl --write   (writes CAPABILITIES.md)
#
# Run once per MAJOR/MINOR release (CLAUDE.md Rule 17), never for a patch release, and never hand-edited
# in between -- the whole point is that CAPABILITIES.md is always reproducible from these two scripts' own
# inputs, exactly like docs/checkCoverage.jl and tools/regenerateRouteAccuracy.jl.
#
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "capabilityNotes.jl"))

struct Entry
    code        ::String
    description ::String
end


# ===== 1. Parse examples/examples.jl into (letter, title, [Entry]) groups =====================================

function  parseExamplesIndex(path::String)

    headerTitles = Dict{String,String}()
    entries      = Entry[]
    lastIdx      = 0

    for  line in readlines(path)
        mHeader = match(r"println\(\"([A-Za-z])\)\s+(.*?)\"\)", line)
        if  mHeader !== nothing
            headerTitles[mHeader.captures[1]] = strip(mHeader.captures[2])
            lastIdx = 0
            continue
        end
        mEntry = match(r"println\(\"([A-Za-z]{2,3})\)\s?(.*?)\"\)", line)
        if  mEntry !== nothing
            push!(entries, Entry(mEntry.captures[1], strip(mEntry.captures[2])))
            lastIdx = length(entries)
            continue
        end
        mCont = match(r"println\(\"(\s+.*?)\"\)", line)
        if  mCont !== nothing && lastIdx > 0
            e = entries[lastIdx]
            entries[lastIdx] = Entry(e.code, e.description * " " * strip(mCont.captures[1]))
            continue
        end
        lastIdx = 0    # any other line (separator, comment, blank) ends the current continuation run
    end

    byLetter = Dict{String, Vector{Entry}}()
    for  e in entries
        gl = string(e.code[1])
        push!(get!(byLetter, gl, Entry[]), e)
    end

    titleOf(gl) = get(headerTitles, gl, get(GROUP_TITLES, gl, "(untitled in examples.jl)"))

    letters = sort(collect(keys(byLetter)))

    return( [ (letter = gl, title = titleOf(gl), entries = byLetter[gl])  for gl in letters ] )
end


# ===== 2. For one branch code, find the best (most recent, genuinely parseable) `Last successful` date ========

function  bestDate(code::String)

    path = joinpath(ROOT, "examples", "example-$code.jl")
    isfile(path) || return( (status = :nofile, date = nothing) )

    dates = Date[]
    for  line in readlines(path)
        m = match(r"#\s*Last successful:\s*(.*)", line)
        m === nothing && continue
        tail = m.captures[1]
        dm = match(r"(\d{1,2})[-\s]?([A-Za-z]{3})[-\s]?(\d{4})", tail)
        dm === nothing && continue
        try
            push!(dates, Date("$(dm.captures[1])-$(dm.captures[2])-$(dm.captures[3])", dateformat"d-u-y"))
        catch
            # an unparseable date token (shouldn't happen given the regex) is simply not counted
        end
    end

    isempty(dates) && return( (status = :undated, date = nothing) )

    return( (status = :dated, date = maximum(dates)) )

end


# ===== 3. Module-name extraction from the description text, overrides applied last =============================

function  extractModule(code::String, description::String)

    haskey(MODULE_OVERRIDE, code) && return( MODULE_OVERRIDE[code] )

    m = match(r"\bthe\s+([A-Z][A-Za-z0-9]*(?:\.[A-Z][A-Za-z0-9]*)?)\s+module\b", description)
    m !== nothing && return( m.captures[1] )

    m = match(r"\(module\s+([A-Z][A-Za-z0-9]*)\)", description)
    m !== nothing && return( m.captures[1] )

    m = match(r"\b(Cascade|Plasma|Empirical)\.([A-Za-z0-9]+)", description)
    m !== nothing && return( "$(m.captures[1]).$(m.captures[2])" )

    return( nothing )

end


# ===== 4. Render one entry, one group, and the whole file ======================================================

function  renderEntry(io::IO, e::Entry)

    retired = occursin(r"^RETIRED"i, e.description)
    desc    = e.description

    println(io, "### $(e.code)) $(desc)")
    println(io)

    if  retired
        println(io, "- **RETIRED** -- kept here only so a search does not come up empty; follow the pointer above.")
        println(io)
        return( nothing )
    end

    modName = extractModule(e.code, desc)
    modName !== nothing  &&  println(io, "- Module: `$modName`")

    d = bestDate(e.code)
    if      d.status == :nofile
        println(io, "- No dedicated `examples/example-$(e.code).jl` file was found for this branch code.")
    elseif  d.status == :undated
        println(io, "- Start here: `examples/example-$(e.code).jl`")
        println(io, "- Date: **unknown/undated** -- the branch exists and runs; nobody has yet recorded a judged date.")
    else
        println(io, "- Start here: `examples/example-$(e.code).jl`, dated **$(Dates.format(d.date, "dd-u-yyyy"))**")
    end

    kw = get(SYNONYMS, e.code, String[])
    !isempty(kw)  &&  println(io, "- Keywords: " * join(kw, ", "))

    restr = get(RESTRICTIONS, e.code, "")
    !isempty(restr)  &&  println(io, "- Restriction: $restr")

    println(io)

    return( nothing )

end

function  renderGroup(io::IO, g)

    println(io, "## $(g.letter)) $(g.title)")
    println(io)
    for  e in g.entries
        renderEntry(io, e)
    end

    return( nothing )

end

function  buildRegister()

    groups = parseExamplesIndex(joinpath(ROOT, "examples", "examples.jl"))

    io = IOBuffer()
    println(io, "# JAC capabilities register")
    println(io)
    println(io, "**What this file is for.** JAC has grown past the point where any one person -- maintainer, " *
                "collaborator, or an AI tool reading the repository -- can hold its full feature set in their " *
                "head. This file exists to answer two questions quickly, for a human or for an AI assistant " *
                "acting on a human's behalf:")
    println(io)
    println(io, "1. **Is JAC already able to do what I need, even approximately?**")
    println(io, "2. **Which one or two example branches should I open first to try it and adapt it?**")
    println(io)
    println(io, "**What this file is NOT.** It is not the manual, and completeness of DETAIL is explicitly not " *
                "the goal -- simple use and fast orientation are. The full technical detail lives in each " *
                "module's own docstrings and in `docs/`; this file only tells you whether to look there at all, " *
                "and where to start. It IS complete in COVERAGE, though: every branch listed by " *
                "`examples/examples.jl` gets its own row here, mechanically, so nothing in that index goes " *
                "missing from this one.")
    println(io)
    println(io, "**How to trust what it says.** Every row's date comes from that branch's own " *
                "`# Last successful: DD-Mon-YYYY` comment line in `examples/*.jl` -- the one marker in this " *
                "codebase that a human, not a tool, has ever written after physically judging the numbers " *
                "(see `CLAUDE.md`, Rule 7). A row marked **unknown/undated** means exactly that: the module " *
                "exists and the branch runs, but nobody has yet judged the physics, or an earlier judgement " *
                "was withdrawn. An omitted `Module:` line means the same kind of honesty -- this register's " *
                "generator could not confidently name one, rather than guessing.")
    println(io)
    println(io, "**Restrictions column.** A one-line hint, not a specification. It exists to stop a false " *
                "start, not to replace reading the branch. A blank entry means no restriction is currently " *
                "flagged -- not that none exists.")
    println(io)
    println(io, "**One file, several branches.** Several example files carry more than one `if true`/`if " *
                "false` branch, each dated independently. The date shown below is the MOST RECENT judged " *
                "date found ANYWHERE in that file, which may belong to a different branch than the one you " *
                "actually want -- always check the specific branch's own `Last successful` comment before " *
                "trusting a copied date. `example-Ar.jl` is a real case of this: one branch there carries a " *
                "genuine 13-Sep-2026 date while a different, still-undated branch in the SAME file is the " *
                "one most readers would reach for first.")
    println(io)
    println(io, "**Searching this file.** Use your editor's find, a plain `grep -i -B1 -A5 \"<term>\" " *
                "CAPABILITIES.md`, or the small helper that prints a matched entry with its full context: " *
                "`julia tools/searchCapabilities.jl \"<term>\"`. All three are plain substring search -- " *
                "\"DR\" will also match \"quaDRupole\"; that is a known, accepted limit of keeping the tool " *
                "genuinely simple, not a defect.")
    println(io)
    println(io, "**Once you have found your branch, the next question is how to turn it into a number worth " *
                "quoting** -- and that is a different document. `docs/UserSteps-JAC.txt` walks an occasional " *
                "user through ten steps from a physics question to a publishable result: classifying the " *
                "problem against JAC's seven kinds of computation, choosing the mean field and the " *
                "correlation model as a PHYSICS decision rather than a default, running the smallest case " *
                "first and then scaling ONE axis at a time, reading JAC's own printed verdicts before " *
                "reading the number itself, and finally comparing against something JAC did not compute. " *
                "This register answers WHETHER and WHERE; those ten steps answer HOW -- and their step 2 is " *
                "the one that sends a reader back here.")
    println(io)
    println(io, "---")
    println(io)

    for  g in groups
        renderGroup(io, g)
        println(io, "---")
        println(io)
    end

    commit = try
        strip(read(`git -C $ROOT rev-parse --short HEAD`, String))
    catch
        "unknown"
    end
    println(io, "*Mechanically generated by `tools/buildCapabilityRegister.jl` on " *
                "$(Dates.format(Dates.now(), "dd-u-yyyy")) from commit `$commit`. Regenerated once per " *
                "MAJOR/MINOR release (`CLAUDE.md` Rule 17) -- never hand-edited in between.*")

    return( String(take!(io)) )

end

if  abspath(PROGRAM_FILE) == @__FILE__

    text = buildRegister()

    if  "--write" in ARGS
        outPath = joinpath(ROOT, "CAPABILITIES.md")
        write(outPath, text)
        println("Wrote $outPath  ($(count(==('\n'), text)) lines).")
    else
        print(text)
    end

end
