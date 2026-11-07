# Functions in this file cover: the external-format writers and the spin-angular term tables --
# the modules that nothing else in src/ calls, so that nothing else can exercise them.
# Alphabetical order within this file.


"""
`TestFrames.exportTestDataSet()`
    ... builds the one small, hand-made `ExportDataSet.DataSet` that every writer test in this file
        consumes:  three levels of a fictitious C^2+ ion and two transitions between them.  It is
        deliberately NOT taken from a computation, because the writers are being tested and not the
        physics, and a hand-made set lets each field be checked against a value written here.

        ONE TRANSITION IS DELIBERATELY DOWNHILL IN ENERGY (upper index 1, lower index 3), so that
        `ExportChianti.writeWgfa`'s own `deltaE <= 0 && continue` guard is exercised rather than
        merely present:  a writer that dropped that guard would emit a negative wavelength.

        A `data::ExportDataSet.DataSet` is returned.
"""
function exportTestDataSet()
    levels = [ ExportDataSet.LevelRecord(1, "2s^2 1S",   1, 0, 0.0,      0.0),
               ExportDataSet.LevelRecord(2, "2s2p 3P0",  3, 1, 0.0,  52367.1),
               ExportDataSet.LevelRecord(3, "2s2p 1P1",  1, 1, 1.0, 102352.0) ]
    trans  = [ ExportDataSet.TransitionRecord(3, 1, 1.769e+09, [1.23e-01, 4.56e-01, 7.89e-01]),
               ExportDataSet.TransitionRecord(1, 3, 0.0,       [1.00e-03, 2.00e-03, 3.00e-03]) ]

    return( ExportDataSet.DataSet("C", 2, 6, 386241.0, levels, [1.0e4, 1.0e5, 1.0e6], trans) )
end


"""
`TestFrames.testModule_ExportAdas(; short::Bool=true)`
    ... tests on module ExportAdas:  the seven-character adf04 number format and the line structure of
        a written adf04 file.  Nothing in `src/` calls this module, so the suite is the only thing that
        can execute it at all.  A `success::Bool` is returned.
"""
function testModule_ExportAdas(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ExportAdas-new.sum")
    printstyled("\n\nTest the module  ExportAdas  ... \n", color=:cyan)
    success = true
    data    = exportTestDataSet()

    # THE adf04 NUMBER FORMAT IS EXACTLY SEVEN CHARACTERS, mantissa + sign + two-digit exponent, and
    # the format's own column layout depends on it;  a number one character wide or narrow shifts every
    # field after it on the line.  Its own docstring states the seven, so this checks the stated contract.
    # THE DOMAIN IS THE NON-NEGATIVE ONE, and the function now RAISES outside it -- the maintainer's
    # decision of 07-Nov-2026, after this very test found that a negative argument returned EIGHT
    # characters ('-3.21+04') against the seven the docstring promised.  The seven-character adf04
    # field has no column for a sign, so a negative would shift every field after it on the line.
    for  x  in  [0.0, 1.0, 1.234e-05, 9.87e+12, 3.21e+04, 1.0e-99, 5.0e+99]
        sa = ExportAdas.formatExponential(x)
        if  length(sa) != 7     success = false
            println("** ExportAdas.formatExponential($x) = '$sa' is $(length(sa)) characters, not 7")   end
    end
    if  ExportAdas.formatExponential(0.0) != "0.00+00"
        success = false;   println("** ExportAdas.formatExponential(0.0) is not '0.00+00'")             end
    # AND THE GUARD MUST FIRE.  A test that only checks the non-negative domain would pass just as
    # well against the old, silently-misaligning version, so this is the assertion that pins the fix.
    try    ExportAdas.formatExponential(-3.21e+04)
           success = false;   println("** formatExponential accepted a negative argument")
    catch  ex
        if  !isa(ex, ErrorException)
            success = false;   println("** formatExponential(-3.21e+04) raised $(typeof(ex)), not an error")  end
    end

    # THE MANTISSA-CARRY BRANCH:  9.999e+03 rounds to 10.00, which must become 1.00+04 and not 10.00+03.
    if  ExportAdas.formatExponential(9.999e+03) != "1.00+04"
        success = false
        println("** the mantissa carry is wrong: 9.999e+03 gives '$(ExportAdas.formatExponential(9.999e+03))'")  end

    # THE FILE STRUCTURE:  one header line, one line per level, a -1 terminator, one temperature line,
    # one line per transition, a second -1.  A writer that lost a terminator would still produce output.
    io = IOBuffer();   ExportAdas.writeAdf04(io, data)
    lines = split(strip(String(take!(io))), "\n")
    expect = 1 + length(data.levels) + 1 + 1 + length(data.transitions) + 1
    if  length(lines) != expect
        success = false;   println("** writeAdf04 wrote $(length(lines)) lines, expected $expect")      end
    if  count(l -> strip(l) == "-1", lines) != 2
        success = false;   println("** writeAdf04 did not write exactly two -1 terminators")            end
    if  !occursin("C++", lines[1])
        success = false;   println("** writeAdf04's header does not name the ion: '$(lines[1])'")        end

    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Format and structure checks only for ExportAdas; NO physics is compared against approved data.")
    testPrint("testModule_ExportAdas()::", success)
    return( success )
end


"""
`TestFrames.testModule_ExportChianti(; short::Bool=true)`
    ... tests on module ExportChianti:  the line structure of a written `.elvlc` file, and the oscillator
        strength that `.wgfa` derives from the A-value, which is re-derived here independently of the
        module.  A `success::Bool` is returned.
"""
function testModule_ExportChianti(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ExportChianti-new.sum")
    printstyled("\n\nTest the module  ExportChianti  ... \n", color=:cyan)
    success = true
    data    = exportTestDataSet()

    io = IOBuffer();   ExportChianti.writeElvlc(io, data)
    lines = split(strip(String(take!(io))), "\n")
    if  length(lines) != length(data.levels) + 1
        success = false;   println("** writeElvlc wrote $(length(lines)) lines for $(length(data.levels)) levels")  end
    if  strip(lines[end]) != "-1"
        success = false;   println("** writeElvlc did not terminate with -1")                            end

    # ONLY THE UPHILL TRANSITION MAY APPEAR.  The test data carries one transition with deltaE <= 0 on
    # purpose;  emitting it would give a negative wavelength, which is why the writer skips it.
    io = IOBuffer();   ExportChianti.writeWgfa(io, data)
    lines = split(strip(String(take!(io))), "\n")
    uphill = count(tr -> data.levels[tr.upperIndex].energy > data.levels[tr.lowerIndex].energy, data.transitions)
    if  length(lines) != uphill + 1
        success = false
        println("** writeWgfa wrote $(length(lines)) lines for $uphill uphill transition(s) plus a terminator")  end

    # THE gf CONVERSION, re-derived here rather than read from the module:  gf = g_upper A lambda^2 / 6.670e15
    # with lambda in Angstrom.  This is the one number in the file that is computed rather than copied.
    if  length(lines) >= 2
        tr     = data.transitions[1]
        upper  = data.levels[tr.upperIndex];     lower = data.levels[tr.lowerIndex]
        lambda = 1.0e8 / (upper.energy - lower.energy)
        gfRef  = (2*upper.J + 1) * tr.aValue * lambda^2 / 6.670e15
        tokens = split(lines[1])
        gfFile = length(tokens) >= 4 ? parse(Float64, tokens[4]) : 0.0
        if  abs(gfFile - gfRef) > 1.0e-4 * abs(gfRef)
            success = false;   println("** writeWgfa's gf = $gfFile against an independent $gfRef")      end
    end

    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Format and structure checks for ExportChianti, plus an independent gf re-derivation.")
    testPrint("testModule_ExportChianti()::", success)
    return( success )
end


"""
`TestFrames.testModule_ExportDataSet(; short::Bool=true)`
    ... tests on module ExportDataSet:  the three container types construct, carry their fields unchanged,
        and accept the empty case that a computation yielding no transitions would hand them.
        A `success::Bool` is returned.
"""
function testModule_ExportDataSet(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ExportDataSet-new.sum")
    printstyled("\n\nTest the module  ExportDataSet  ... \n", color=:cyan)
    success = true
    data    = exportTestDataSet()

    if  data.element != "C"  ||  data.ionCharge != 2  ||  data.nuclearCharge != 6
        success = false;   println("** ExportDataSet.DataSet did not carry its own identification")     end
    if  length(data.levels) != 3  ||  length(data.transitions) != 2
        success = false;   println("** ExportDataSet.DataSet did not carry its own records")            end
    if  data.levels[3].label != "2s2p 1P1"  ||  data.levels[3].energy != 102352.0
        success = false;   println("** ExportDataSet.LevelRecord did not carry its own fields")         end
    # EVERY upsilon LIST MUST ALIGN ONE-TO-ONE WITH THE OWNING SET'S TEMPERATURES -- the docstring says so,
    # and a writer indexes one by the other, so a mismatch here is an out-of-bounds error at write time.
    for  tr  in  data.transitions
        if  length(tr.upsilons) != length(data.temperatures)
            success = false
            println("** $(length(tr.upsilons)) upsilons against $(length(data.temperatures)) temperatures")  end
    end
    # the empty case: a computation that found no transitions must still give a writable set
    empty = ExportDataSet.DataSet("H", 0, 1, 109678.8, ExportDataSet.LevelRecord[], Float64[],
                                  ExportDataSet.TransitionRecord[])
    io = IOBuffer()
    try    ExportAdas.writeAdf04(io, empty)
    catch  ex    success = false;   println("** an empty DataSet cannot be written: $ex")               end

    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Container checks only for ExportDataSet; it carries no physics by design.")
    testPrint("testModule_ExportDataSet()::", success)
    return( success )
end


"""
`TestFrames.testModule_ExportVamdc(; short::Bool=true)`
    ... tests on module ExportVamdc:  every element the minimal XSAMS writer opens is closed again, and
        one state identifier per level is emitted.  The document is NOT schema-valid by the module's own
        admission, so nothing here checks it against XSAMS.  A `success::Bool` is returned.
"""
function testModule_ExportVamdc(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-ExportVamdc-new.sum")
    printstyled("\n\nTest the module  ExportVamdc  ... \n", color=:cyan)
    success = true
    data    = exportTestDataSet()

    io = IOBuffer();   ExportVamdc.writeMinimalXsams(io, data)
    xml = String(take!(io))

    # TAG BALANCE IS THE ONE THING WORTH CHECKING without a schema: an unclosed element makes the file
    # unreadable by any XML parser, and is exactly what a hand-written writer gets wrong when a loop
    # gains an early `continue`.
    for  tag  in  ["XSAMSData", "Species", "Atom", "AtomicStates", "RadiativeTransitions"]
        nOpen  = count(i -> true, eachmatch(Regex("<" * tag * "[ >]"), xml))
        nClose = count(i -> true, eachmatch(Regex("</" * tag * ">"),   xml))
        if  nOpen != nClose  ||  nOpen == 0
            success = false;   println("** <$tag> opened $nOpen times and closed $nClose times")        end
    end
    nStates = count(i -> true, eachmatch(r"<AtomicState stateID=", xml))
    if  nStates != length(data.levels)
        success = false;   println("** $nStates AtomicState elements for $(length(data.levels)) levels") end
    nTrans  = count(i -> true, eachmatch(r"<RadiativeTransition>", xml))
    if  nTrans != length(data.transitions)
        success = false
        println("** $nTrans RadiativeTransition elements for $(length(data.transitions)) transitions")  end

    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Tag-balance checks only for ExportVamdc; the document is not schema-valid by design.")
    testPrint("testModule_ExportVamdc()::", success)
    return( success )
end


"""
`TestFrames.testModule_SpinAngularTables(; short::Bool=true)`
    ... tests on module SpinAngularTables:  the quasispin projection and the quasispin itself, each against
        a law that holds whatever sign convention the module uses, so that the test survives a convention
        change and fails on a defect.  A `success::Bool` is returned.
"""
function testModule_SpinAngularTables(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-SpinAngularTables-new.sum")
    printstyled("\n\nTest the module  SpinAngularTables  ... \n", color=:cyan)
    success = true

    # The module is mentioned in ONE comment in src/ and called from nowhere, so the suite is the only
    # thing that executes it.  Two laws are checked, both convention-free:
    #
    #   the quasispin PROJECTION M_Q = (N - Omega)/2 is antisymmetric about half filling, so that
    #   M_Q(N) = -M_Q(2 Omega - N) and M_Q = 0 at N = Omega;  and
    #
    #   the quasispin Q = (Omega - nu)/2 falls by 1/2 per unit of seniority and reaches 0 at nu = Omega,
    #   where Omega = (2j+1)/2 is half the subshell's capacity.
    #
    # Neither law names a sign convention, so a change of convention keeps them while a wrong index, a
    # wrong capacity or an off-by-one in the seniority breaks them.
    for  j  in  [AngularJ64(1//2), AngularJ64(3//2), AngularJ64(5//2), AngularJ64(7//2), AngularJ64(9//2)]
        capacity = Basics.twice(j) + 1
        omega    = div(capacity + 1, 2)
        for  N  in  0:capacity
            mLeft  = SpinAngularTables.qshellTermM(j, N)
            mRight = SpinAngularTables.qshellTermM(j, capacity - N)
            if  Basics.twice(mLeft) != -Basics.twice(mRight)
                success = false
                println("** M_Q is not antisymmetric for j = $j, N = $N:  $mLeft against $mRight")      end
        end
        if  Basics.twice(SpinAngularTables.qshellTermM(j, omega)) != 0
            success = false;   println("** M_Q is not zero at half filling for j = $j")                 end
        qPrev = nothing
        for  nu  in  0:omega
            q = SpinAngularTables.qshellTermQ(j, nu)
            if  !isnothing(qPrev)  &&  Basics.twice(qPrev) - Basics.twice(q) != 1
                success = false
                println("** Q does not fall by 1/2 per seniority for j = $j, nu = $nu")                 end
            qPrev = q
        end
        if  Basics.twice(SpinAngularTables.qshellTermQ(j, omega)) != 0
            success = false;   println("** Q is not zero at nu = Omega for j = $j")                      end
        # THE TWO FUNCTIONS MUST AGREE WITH EACH OTHER at the closed shell: Q(nu = 0) = -M_Q(N = 0),
        # since both equal Omega/2 in magnitude.  This is what catches a capacity used in one and not
        # the other, which neither law above would see on its own.
        if  Basics.twice(SpinAngularTables.qshellTermQ(j, 0)) != -Basics.twice(SpinAngularTables.qshellTermM(j, 0))
            success = false;   println("** Q(nu=0) and M_Q(N=0) disagree in magnitude for j = $j")       end
    end

    Defaults.setDefaults("print summary: close", "")
    printTest, iostream = Defaults.getDefaults("test flag/stream")
    println(iostream, "Quasispin laws for SpinAngularTables; no coefficient of fractional parentage is compared.")
    testPrint("testModule_SpinAngularTables()::", success)
    return( success )
end
