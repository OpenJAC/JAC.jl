
"""
`Basics.display()`  
    ... reports/displays various information about the (present) settings or some given data, often in some separate window; 
        nothing is returned if not indicated otherwise. Cf. JAC.Defaults.setDefaults(). 

  + `("constants")`  or  `("physical constants")`  ... to display (all) currently defined physical constants.

  + `("settings")`     ... to display (all) currently defined settings of the JAC module.
"""
function Basics.display(sa::String)

    if        sa in ["constants", "physical constants"]
        println("Physical constants are defines as follows:  \n", 
                "------------------------------------------  \n")
        sb = "  + Fine-structure constant:    " * string( Defaults.getDefaults("alpha") );                println(sb)        
        sb = "  + Electron mass [kg]:         " * string( Defaults.getDefaults("electron mass: kg") );    println(sb)        
        sb = "  + Electron mass [amu]:        " * string( Defaults.getDefaults("electron mass: amu") );   println(sb)        
        println()

    elseif   sa == "settings"
        println("Current settings of the JAC module:  \n", 
                "-----------------------------------  \n")
        sb = "  + Framework:                              " * string( Defaults.getDefaults("framework") );                println(sb)        
        sb = "  + Energy unit:                            " * string( Defaults.getDefaults("unit: energy") );             println(sb)        
        sb = "  + Rate and transition probability unit:   " * string( Defaults.getDefaults("unit: rate") );               println(sb)        
        sb = "  + Cross section unit:                     " * string( Defaults.getDefaults("unit: cross section") );      println(sb) 
        sb = "  + Time unit:                              " * string( Defaults.getDefaults("unit: time") );               println(sb) 
        println()       
        
        if      Defaults.getDefaults("standard grid") != false    println("  + A standard grid has been defined; cf. JAC.Defaults.getDefaults()" )
        elseif  Defaults.getDefaults("standard grid") == false    println("  + No standard grid has yet been defined; cf. JAC.Defaults.setDefaults()" )
        end
        println()       
 
    else    error("Unsupported keystring:: $sa")
    end

    return( nothing )
end
