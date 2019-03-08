

# export  warn

"""
`JAC.warn()`  ... deals with warnings that occur during a run and session; it handles the global array JAC_WARNINGS.

  + `(AddWarning, warning::String)`  ... to add warning to the global array JAC_WARNINGS.

  + `(PrintWarnings)`  ... to print all warnings that are currently kept in the global array JAC_WARNINGS.

  + `(ResetWarnings)`  ... to reset the global array JAC_WARNINGS.
"""
function warn(wa::Warnings)
    global JAC_WARNINGS

    if        wa == PrintWarnings
        iostream = open("jac-warn.report", "w") 
        println(iostream, " ")
        println(iostream, "\n\n",
                          "Warnings from the present sessions or run, printed at $( string(now())[1:16] ): \n",
                          "======================================================================= \n")
        for sa in  JAC_WARNINGS
            println(iostream, "++ " * sa)
        end
        close(iostream)
        #
    elseif    wa == ResetWarnings
        ## @warn("Reset global array JAC_WARNINGS.")
        printstyled("JAC.warn():  Reset global array JAC_WARNINGS.", color=:light_magenta)
        JAC_WARNINGS = String[]
    else      error("Unsupported Warnings:: $wa")
    end

    nothing
end


function warn(wa::Warnings, sa::String)
    global JAC_WARNINGS

    if        wa == AddWarning
        push!(JAC_WARNINGS, sa)
    else      error("Unsupported Warnings:: $wa")
    end

    nothing
end
