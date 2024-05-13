
println("Cj) Apply & test the MultipolePolarizibility module with ASF from an internally generated multiplet.")

if  true
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-MultipolePolarizibility.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            properties=[Polarity()], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            polaritySettings=MultipolePolarizibility.Settings(EmMultipole[], 0, 0, Float64[], false, false, Int64[]) )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end


