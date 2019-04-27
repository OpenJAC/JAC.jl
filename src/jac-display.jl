
## export  display

"""
`JAC.Basics.display()`  ... reports/displays various information about the (present) settings or some given data, often in some separate window; 
                     nothing is returned if not indicated otherwise. Cf. JAC.Constants.define(). 

  + `("constants")`  or  `("physical constants")`  ... to display (all) currently defined physical constants.

  + `("settings")`     ... to display (all) currently defined settings of the JAC module.
"""
function display(sa::String)

    if        sa in ["constants", "physical constants"]
        println("Physical constants are defines as follows:  \n", 
                "------------------------------------------  \n")
        sb = "  + Fine-structure constant:    " * string( Constants.give("alpha") );                println(sb)        
        sb = "  + Electron mass [kg]:         " * string( Constants.give("electron mass: kg") );    println(sb)        
        sb = "  + Electron mass [amu]:        " * string( Constants.give("electron mass: amu") );   println(sb)        
        println()

    elseif   sa == "settings"
        println("Current settings of the JAC module:  \n", 
                "-----------------------------------  \n")
        sb = "  + Framework:                              " * string( Constants.give("framework") );                println(sb)        
        sb = "  + Energy unit:                            " * string( Constants.give("unit: energy") );             println(sb)        
        sb = "  + Rate and transition probability unit:   " * string( Constants.give("unit: rate") );               println(sb)        
        sb = "  + Cross section unit:                     " * string( Constants.give("unit: cross section") );      println(sb) 
        sb = "  + Time unit:                              " * string( Constants.give("unit: time") );               println(sb) 
        println()       
        
        if      Constants.give("standard grid") != false    println("  + A standard grid has been defined; cf. JAC.Constants.give()" )
        elseif  Constants.give("standard grid") == false    println("  + No standard grid has yet been defined; cf. JAC.Constants.define()" )
        end
        println()       
 
    else    error("Unsupported keystring:: $sa")
    end

    return( nothing )
end
