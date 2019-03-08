
## export  display

"""
`JAC.display()`  ... reports/displays various information about the (present) settings or some given data, often in some separate window; 
                     nothing is returned if not indicated otherwise. Cf. JAC.define(). 

  + `("constants")`  or  `("physical constants")`  ... to display (all) currently defined physical constants.

  + `("settings")`     ... to display (all) currently defined settings of the JAC module.
"""
function display(sa::String)

    if        sa in ["constants", "physical constants"]
        println("Physical constants are defines as follows:  \n", 
                "------------------------------------------  \n")
        sb = "  + Fine-structure constant:    " * string( JAC.give("alpha") );                println(sb)        
        sb = "  + Electron mass [kg]:         " * string( JAC.give("electron mass: kg") );    println(sb)        
        sb = "  + Electron mass [amu]:        " * string( JAC.give("electron mass: amu") );   println(sb)        
        println()

    elseif   sa == "settings"
        println("Current settings of the JAC module:  \n", 
                "-----------------------------------  \n")
        sb = "  + Framework:                              " * string( JAC.give("framework") );                println(sb)        
        sb = "  + Energy unit:                            " * string( JAC.give("unit: energy") );             println(sb)        
        sb = "  + Rate and transition probability unit:   " * string( JAC.give("unit: rate") );               println(sb)        
        sb = "  + Cross section unit:                     " * string( JAC.give("unit: cross section") );      println(sb) 
        sb = "  + Time unit:                              " * string( JAC.give("unit: time") );               println(sb) 
        println()       
        
        if      JAC.give("standard grid") != false    println("  + A standard grid has been defined; cf. JAC.give()" )
        elseif  JAC.give("standard grid") == false    println("  + No standard grid has yet been defined; cf. JAC.define()" )
        end
        println()       
 
    else    error("Unsupported keystring:: $sa")
    end

    return( nothing )
end
