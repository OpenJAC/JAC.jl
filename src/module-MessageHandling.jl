
"""
`module MessageHandling`  ... a submodel of JAC that contains various methods to deal with messages and warnings. Such messages refer 
                              to warnings, time settings or just limitations of the code and are typically issued in course of a computation 
                              and can later be printed together at well-defined times. 
"""
module MessageHandling

    # using JAC
    messageList    = String[]
    fromMethodList = String[]
    keystringList  = String[]


    """
    `JAC.MessageHandling.empty()`  ... to reset all lists for collecting messages; nothing is returned.
    """
    function  empty()
        messageList = String[];    fromMethodList = String[];    keystringList = String[]
        return( nothing )
    end


    """
    `JAC.MessageHandling.addMessage(fromMethod::String, keystring::String, message::String)`  ... to add one 'message' to the list with 
                                    some additonal information from where it was issued.
    """
    function  addMessage(fromMethod::String, keystring::String, message::String)
        push!(messageList, message);    push!(fromMethodList, fromMethod);    push!(keystringList, keystring);
        return( nothing )
    end


    """
    `JAC.MessageHandling.printMessages()`  ... to print out all current messages in a neat format.
    """
    function printMessages()
        date = string( Base.Dates.now() )
        println("\nMessages on $(date[1:10]) at $(date[12:19]):") 
        println("-----------------------------------") 
        for  i = 1:length(messageList)
            sa = "  + " * messageList[i] * repeat(" ", 100);    sa = sa[1:120] * "[from: " * fromMethodList[i] * ": " * keystringList[i] * "]"
            lprint(sa)
        end

        return( nothing )
    end

end # module
