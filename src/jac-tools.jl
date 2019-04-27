
"""
`Basics.tools()`  ... select different tools from a menu for which no further input is required.
"""
function Basics.tools(x::Int64)
    println("JAC.Tools.select(): A menu if nothing comes in")
    #
    t1 = "Simple tools JAC without input: "
    b1 = dropdown(["Start task", "Grid calculator", "Another task"])
    b2 = slider(1:100, label = "To what extent?", value = 33)
    update = button("Update")
    ui = vbox( hbox( pad(1em, t1) ),
               hbox( pad(1em, b1), pad(1em, b2), pad(1em, update) )
             )
    Interact.display(ui)
    output = Interact.@map (&update;  (observe(b1)[], observe(b2)[]) ) 
    return( output )
end


"""
`Basics.tools(dict::Dict)`  ... select different tools from a menu if proper results::dict are given in terms of a dictionary as obtained
                             usually from an Atomic.Computation.
"""
function Basics.tools(dict::Dict)
    println("JAC.Tools.select(): A menu if a dict comes in; $dict ")
end

