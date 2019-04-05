
## export  test

function testCompareFiles(fold::String, fnew::String, sa::String, noLines::Int64) 
    success = true
    printTest, iostream = JAC.give("test flag/stream")
    
    # Compare the test computations with previous results
    oldLines = readlines(fold)
    newLines = readlines(fnew)
    iold = 0;  for i=1:length(oldLines)   line = oldLines[i];   if  occursin(sa, line)  iold = i;   break   end   end  
    inew = 0;  for i=1:length(newLines)   line = newLines[i];   if  occursin(sa, line)  inew = i;   break   end   end  
    if  iold == 0   ||   inew == 0    success = false 
        if printTest   println(iostream, "Tries to compare two inappropriate files fold = $(fold); fnew =$(fnew) on string $sa  ($iold, $inew)")  end
        return( success )
    end
    #
    ii = inew + 1
    for  i = iold+2:iold+noLines
       ii = ii + 1
       if   length(oldLines[i]) < 5    continue    end
       if   oldLines[i] != newLines[ii]    success = false
           if  printTest  println(iostream, "    *** Old::  " * oldLines[i]) 
                          println(iostream, "    *** New::  " * newLines[ii])
           end
       end
    end
    
    return( success )
end 


function testPrint(sa::String, success::Bool) 
    printTest, iostream = JAC.give("test flag/stream")
    ok(succ) =  succ ? "[OK]" : "[Fail]"
    sb = sa * JAC.TableStrings.hBlank(110);   sb = sb[1:100] * ok(success);    println(iostream, sb)
    return( nothing )
end 



"""
`JAC.testMethod_integrate_ongrid(; short::Bool=true)`  ... tests the integration on grid.
"""
function testMethod_integrate_ongrid(; short::Bool=true)
    success = true
    printTest, iostream = JAC.give("test flag/stream")

    # Test the integration on the grid for an analytical function
    grid = JAC.Radial.Grid("grid: by given parameters", rnt=2.0e-10, h=5.0e-3, NoPoints=9000)
    function f1(r)    
        A = 10^5;    gamma = 3;    a = 300
        return( A * r^gamma * exp(-a * r^2) )
    end
    exact1 = 5/9

    integrand = Vector{Float64}(size(grid.r, 1))
    for i = 1:size(grid.r, 1)
        integrand[i] = f1(grid.r[i])
    end
    integrand = integrand .* grid.rp[1:size(integrand, 1)]   # adopt to the form of Grasp92

    integral  = JAC.integrate("function: on radial grid, Newton-Cotes", integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... Newton-Cotes:  I = $integral,  Err = $err")  end
    end

    integral  = JAC.integrate("function: on radial grid, Simpson rule", integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... Simpson rule:  I = $integral,  Err = $err")  end
    end

    integral  = JAC.integrate("function: on radial grid, trapez rule", integrand, grid)
    err = abs(integral - exact1)
    if  abs(err) > 1.0e-12
        success = false
        if printTest   info(iostream, "... trapez rule:  I = $integral,  Err = $err")  end
    end
    
    println(iostream, "Warning(testMethod_integrate_ongrid): test of integration with Grasp orbitals has been set silent.")
    #=
    # Test the integration with orbital functions from Grasp92
    grid = JAC.Radial.Grid("grid: exponential")

    orbitals1 = JAC.readOrbitalFileGrasp92("../test/approved/Ne-0+-scf.exp.out", grid)
    orbitals2 = JAC.readOrbitalFileGrasp92("../test/approved/Ne-1+-scf.exp.out", grid)
 
    for i = 1:size(orbitals1, 1)
        for j = 1:size(orbitals2, 1)
            orb1      = orbitals1[i];   orb2 = orbitals2[j]
            if   orb1.subshell.kappa != orb2.subshell.kappa    break   end

            mtp       = min(size(orb1.P, 1), size(orb2.P, 1))
            integrand = ( orb1.P[1:mtp] .* orb2.P[1:mtp] + orb1.Q[1:mtp] .* orb2.Q[1:mtp] ) .* grid.rp[1:mtp]
      
            integrala = JAC.integrate("function: on radial grid, Newton-Cotes", integrand, grid)^2
            integralb = JAC.integrate("function: on radial grid, Simpson rule", integrand, grid)^2
            integralc = JAC.integrate("function: on radial grid, trapez rule",  integrand, grid)^2
      
            info(iostream, "<$(string(orb1.subshell)) | $(string(orb2.subshell))> = $integrala, $integralb, $integralc")
        end
    end  =#

    testPrint("testMethod_integrate_ongrid()::", success)
    return(success)  
end


"""
`JAC.testMethod_Wigner_3j(; short::Bool=true)`  ... tests on Wigner 3j symbols.
"""
function testMethod_Wigner_3j(; short::Bool=true)
    success = true
    printTest, iostream = JAC.give("test flag/stream")

    # Wigner_3j(1,2,1,0,0,0) = 0.36514837167011074230
    a = c = AngularJ64(1);    b = AngularJ64(2);    ma = mb = mc = AngularM64(0)
    wa = JAC.AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.36514837167011074230) > 1.0e-12
        success = false
        if printTest   info(iostream, "Wigner_3j(1,2,1,0,0,0) = 0.36514837167011074230 ... but obtains value = $wa")   end
    end

    # Wigner_3j(3,6,3,0,0,0) = 0.182482967150452976281
    a = c = AngularJ64(3);    b = AngularJ64(6);    ma = mb = mc = AngularM64(0)
    wa = JAC.AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.18248296715045297628) > 1.0e-12
        success = false
        if printTest   info(iostream, "Wigner_3j(3,6,3,0,0,0) = 0.182482967150452976281 ... but obtains value = $wa")   end
    end

    # Wigner_3j(1/2,1,1/2,1/2,0,-1/2) = 0.40824829046386301637
    a = c = AngularJ64(1//2);    b = AngularJ64(1);    ma =  AngularM64(1//2);   mb =  AngularM64(0);    mc = AngularM64(-1//2)
    wa = JAC.AngularMomentum.Wigner_3j(a,b,c,ma,mb,mc)
    if  abs(wa - 0.40824829046386301637) > 1.0e-12    
        success = false
        if printTest   info(iostream, "Wigner_3j(1/2,1,1/2,1/2,0,-1/2) = 0.40824829046386301637 ... but obtains value = $wa")  end
    end

    testPrint("testMethod_Wigner_3j()::", success)
    return(success)  
end



"""
`JAC.testModule_AlphaVariation(; short::Bool=true)`  ... tests on module JAC.AlphaVariation.
"""
function testModule_AlphaVariation(; short::Bool=true) 
    JAC.define("print summary: open", "test-AlphaVariation-new.sum")
    printstyled("\n\nTest the module  AlphaVariation  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Model(26.); properties=[JAC.AlphaX], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            alphaSettings=AlphaVariation.Settings(true, true, false, Int64[]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-AlphaVariation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-AlphaVariation-new.sum"), "Alpha variation parameters:", 1) 
    testPrint("testModule_AlphaVariation()::", success)
    return(success)  
end



"""
`JAC.testModule_Auger(; short::Bool=true)`  ... tests on module JAC.Auger.
"""
function testModule_Auger(; short::Bool=true) 
    JAC.define("print summary: open", "test-Auger-new.sum")
    printstyled("\n\nTest the module  Auger  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Model(36.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, NoPoints = 600),  
                            initialConfigs=[Configuration("1s^2 2s^2 2p"), Configuration("1s 2s^2 2p^2")],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")], process = JAC.AugerX,
                            processSettings = Auger.Settings(true, true, true, Tuple{Int64,Int64}[(3,1), (4,1), (5,1), (6,1)], 0., 1.0e6, 2, "Coulomb") )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Auger-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Auger-new.sum"), "Auger rates and intr", 17) 
    testPrint("testModule_Auger()::", success)
    return(success)  
end



"""
`JAC.testModule_CoulombExcitation(; short::Bool=true)`  ... tests on module JAC.CoulombExcitation.
"""
function testModule_CoulombExcitation(; short::Bool=true) 
    JAC.define("print summary: open", "test-CoulombExcitation-new.sum")
    printstyled("\n\nTest the module  CoulombExcitation  ... \n", color=:cyan)
    ### Make the tests
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = JAC.give("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-CoulombExcitation-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Auger-approved.sum"), 
    ##                             joinpath(@__DIR__, "..", "test", "test-Auger-new.sum"), "Auger rates and intrinsic angular parameters:", 25) 
    testPrint("testModule_CoulombExcitation()::", success)
    return(success)  
end



"""
`JAC.testModule_DecayYield(; short::Bool=true)`  ... tests on module JAC.DecayYield.
"""
function testModule_DecayYield(; short::Bool=true) 
    JAC.define("print summary: open", "test-DecayYield-new.sum")
    printstyled("\n\nTest the module  DecayYield  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Yields], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            yieldSettings=DecayYield.Settings(true, true, false, Int64[]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-DecayYield-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-DecayYield-new.sum"), "Fluorescence and Auger decay yields:", 6) 
    testPrint("testModule_DecayYield()::", success)
    return(success)  
end



"""
`JAC.testModule_Dielectronic(; short::Bool=true)`  ... tests on module JAC.Dielectronic.
"""
function testModule_Dielectronic(; short::Bool=true) 
    JAC.define("print summary: open", "test-Dielectronic-new.sum")
    printstyled("\n\nTest the module  Dielectronic  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                            process = JAC.Dierec, 
                            processSettings=Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], false, 
                                                                true, Tuple{Int64,Int64,Int64}[(1,1,0)], 0., 0., 0., "Coulomb")  )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Dielectronic-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Dielectronic-new.sum"), 
                               "Partial (Auger) capture", 10) 
    testPrint("testModule_Dielectronic()::", success)
    return(success)  
end



"""
`JAC.testModule_Einstein(; short::Bool=true)`  ... tests on module JAC.Einstein.
"""
function testModule_Einstein(; short::Bool=true) 
    JAC.define("print summary: open", "test-Einstein-new.sum")
    printstyled("\n\nTest the module  Einstein  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(36.); properties=[JAC.EinsteinX], 
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            einsteinSettings=Einstein.Settings([E1, M1, E2, M2], true, 
                            true, Tuple{Int64,Int64}[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0)], 0., 0., 10000. ) )

    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Einstein-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Einstein-new.sum"), "Einstein coefficients, t", 100) 
    testPrint("testModule_Einstein()::", success)
    return(success)  
end



"""
`JAC.testModule_FormFactor(; short::Bool=true)`  ... tests on module JAC.FormFactor.
"""
function testModule_FormFactor(; short::Bool=true) 
    JAC.define("print summary: open", "test-FormFactor-new.sum")
    printstyled("\n\nTest the module  FormFactor  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.FormF], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            formSettings=FormFactor.Settings([0.1], true, false, Int64[]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-FormFactor-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-FormFactor-new.sum"), "Standard and modifi", 6) 
    testPrint("testModule_FormFactor()::", success)
    return(success)  
end



"""
`JAC.testModule_GreenFunction(; short::Bool=true)`  ... tests on module JAC.GreenFunction.
"""
function testModule_GreenFunction(; short::Bool=true) 
    JAC.define("print summary: open", "test-GreenFunction-new.sum")
    printstyled("\n\nTest the module  GreenFunction  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Green], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            greenSettings=GreenFunction.Settings(true, true, false, Int64[]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-GreenFunction-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-GreenFunction-new.sum"), "Green's function summary:", 6) 
    testPrint("testModule_GreenFunction()::", success)
    return(success)  
end



"""
`JAC.testModule_Hfs(; short::Bool=true)`  ... tests on module JAC.Hfs.
"""
function testModule_Hfs(; short::Bool=true) 
    JAC.define("print summary: open", "test-Hfs-b-new.sum")
    printstyled("\n\nTest the module  Hfs  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0); properties=[JAC.HFS],
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            hfsSettings=Hfs.Settings(true, true, true, true, true, true, false, Int64[] ) )

    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Hfs-b-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Hfs-b-new.sum"), "Level        J^P           Energy", 20) 
    testPrint("testModule_Hfs()::", success)
    return(success)  
end



"""
`JAC.testModule_IsotopeShift(; short::Bool=true)`  ... tests on module JAC.IsotopeShift.
"""
function testModule_IsotopeShift(; short::Bool=true) 
    JAC.define("print summary: open", "test-IsotopeShift-new.sum")
    printstyled("\n\nTest the module  IsotopeShift  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Isotope], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            isotopeSettings=IsotopeShift.Settings(true, false, false, true, false, Int64[], "method-1") )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-IsotopeShift-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-IsotopeShift-new.sum"), "IsotopeShift parameters and amplitudes:", 15) 
    testPrint("testModule_IsotopeShift()::", success)
    return(success)  
end



"""
`JAC.testModule_LandeZeeman(; short::Bool=true)`  ... tests on module JAC.LandeZeeman.
"""
function testModule_LandeZeeman(; short::Bool=true) 
    JAC.define("print summary: open", "test-LandeZeeman-new.sum")
    printstyled("\n\nTest the module  LandeZeeman  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26., "Fermi", 58., 3.75, AngularJ64(5//2), 1.0, 2.0); 
                            properties=[JAC.LandeJ], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            zeemanSettings=LandeZeeman.Settings(true, true, true, true, 0., true, false, Int64[] ) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-LandeZeeman-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-LandeZeeman-new.sum"), "Lande g_J factors and Zeeman amplitudes:", 30) 
    testPrint("testModule_LandeZeeman()::", success)
    return(success)  
end



"""
`JAC.testModule_MultiPhotonDeExcitation(; short::Bool=true)`  ... tests on module JAC.MultiPhotonDeExcitation.
"""
function testModule_MultiPhotonDeExcitation(; short::Bool=true) 
    JAC.define("print summary: open", "test-MultiPhotonDeExcitation-new.sum")
    printstyled("\n\nTest the module  MultiPhotonDeExcitation  ... \n", color=:cyan)
    ### Make the tests
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = JAC.give("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-MultiPhotonDeExcitation-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultiPhotonDeExcitation-approved.sum"), 
    ##                             joinpath(@__DIR__, "..", "test", "test-MultiPhotonDeExcitation-new.sum"), "xxx", 100) 
    testPrint("testModule_MultiPhotonDeExcitation()::", success)
    return(success)  
end



"""
`JAC.testModule_MultipoleMoment(; short::Bool=true)`  ... tests on module JAC.MultipoleMoment.
"""
function testModule_MultipoleMoment(; short::Bool=true) 
    JAC.define("print summary: open", "test-MultipoleMoment-new.sum")
    ### Make the tests
    printstyled("\n\nTest the module  MultipoleMoment  ... \n", color=:cyan)
    grid = Radial.Grid("grid: exponential")
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=grid, properties=JAC.AtomicLevelProperty[],
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")] )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    flow = 6;    fup = 8;   ilow = 1;   iup = 3
    println("\n\nDipole amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                MultipoleMoment.dipoleAmplitude(finalLevel, initialLevel, grid; display=true)
            end
        end
    end
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultipoleMoment-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-MultipoleMoment-new.sum"), "Dipole amplitude", 6) 
    testPrint("testModule_MultipoleMoment()::", success)
    return(success)  
end



"""
`JAC.testModule_MultipolePolarizibility(; short::Bool=true)`  ... tests on module JAC.MultipolePolarizibility.
"""
function testModule_MultipolePolarizibility(; short::Bool=true) 
    JAC.define("print summary: open", "test-MultipolePolarizibility-new.sum")
    printstyled("\n\nTest the module  MultipolePolarizibility  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Polarity], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            polaritySettings=MultipolePolarizibility.Settings(EmMultipole[], 0, 0, Float64[], false, false, Int64[]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-MultipolePolarizibility-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-MultipolePolarizibility-new.sum"), 
                               "Multipole polarizibilities and amplitudes:", 5) 
    testPrint("testModule_MultipolePolarizibility()::", success)
    return(success)  
end



"""
`JAC.testModule_ParityNonConservation(; short::Bool=true)`  ... tests on module JAC.ParityNonconservation.
"""
function testModule_ParityNonConservation(; short::Bool=true) 
    JAC.define("print summary: open", "test-ParityNonConservation-new.sum")
    printstyled("\n\nTest the module  ParityNonConservation  ... \n", color=:cyan)
    ### Make the tests
    grid = give("standard grid")
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=JAC.AtomicLevelProperty[],
                            configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")] )

    wxa  = perform(wa; output=true)
    wma  = wxa["multiplet:"]

    nModel = Nuclear.Model(26.);    flow = 6;    fup = 8;   ilow = 1;   iup = 3

    println("\n\nWeak-charge and Schiff-moment amplitudes:\n")
    for  finalLevel in wma.levels
        for  initialLevel in wma.levels
            if  flow <= finalLevel.index <= fup   &&    ilow <= initialLevel.index <= iup
                ParityNonConservation.weakChargeAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
                ParityNonConservation.schiffMomentAmplitude(finalLevel, initialLevel, nModel, grid; display=true)
            end
        end
    end
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-ParityNonConservation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-ParityNonConservation-new.sum"), "weak-charge amplitude", 12) 
    testPrint("testModule_ParityNonConservation()::", success)
    return(success)  
end



"""
`JAC.testModule_PlasmaShift(; short::Bool=true)`  ... tests on module JAC.PlasmaShift.
"""
function testModule_PlasmaShift(; short::Bool=true) 
    JAC.define("print summary: open", "test-PlasmaShift-new.sum")
    printstyled("\n\nTest the module  PlasmaShift  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); properties=[JAC.Plasma], 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            plasmaSettings=PlasmaShift.Settings() )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PlasmaShift-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PlasmaShift-new.sum"), "Plasma screening", 3) 
    testPrint("testModule_PlasmaShift()::", success)
    return(success)  
end



"""
`JAC.testModule_PhotoExcitation(; short::Bool=true)`  ... tests on module JAC.PhotoExcitation.
"""
function testModule_PhotoExcitation(; short::Bool=true) 
    JAC.define("print summary: open", "test-PhotoExcitation-new.sum")
    printstyled("\n\nTest the module  PhotoExcitation  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                            initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")], 
                            process = JAC.PhotoExc, 
                            processSettings=PhotoExcitation.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, true, true, 
                                                                 false, Tuple{Int64,Int64}[], 0., 0., 1.0e6, JAC.ExpStokes(0., 0., 0.) ) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitation-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitation-new.sum"),  
                               "Photo-excitation cross sections for (completely)", 200) 
    testPrint("testModule_PhotoExcitation()::", success)
    return(success)  
end



"""
`JAC.testModule_PhotoExcitationAutoion(; short::Bool=true)`  ... tests on module JAC.PhotoExcitationAutoion.
"""
function testModule_PhotoExcitationAutoion(; short::Bool=true) 
    JAC.define("print summary: open", "test-PhotoExcitationAutoion-new.sum")
    printstyled("\n\nTest the module  PhotoExcitationAutoion  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, NoPoints = 600), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2"), Configuration("1s 2p^2")],
                            finalConfigs=[Configuration("1s^2")], 
                            process = JAC.PhotoExcAuto, 
                            processSettings=PhotoExcitationAutoion.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, Tuple{Int64,Int64,Int64}[(1,1,1)], 2)  )

wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = JAC.give("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-PhotoExcitationAutoion-new.sum")
    success = true
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitationAutoion-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitationAutoion-new.sum"), "xxx", 10) 
    testPrint("testModule_PhotoExcitationAutoion()::", success)
    return(success)  
end



"""
`JAC.testModule_PhotoExcitationFluores(; short::Bool=true)`  ... tests on module JAC.PhotoExcitationFluores.
"""
function testModule_PhotoExcitationFluores(; short::Bool=true) 
    JAC.define("print summary: open", "test-PhotoExcitationFluores-new.sum")
    printstyled("\n\nTest the module  PhotoExcitationFluores  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(26.), 
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                            intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                            finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                            process = JAC.PhotoExcFluor, 
                            processSettings=PhotoExcitationFluores.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], 
                                                                            true, true, Tuple{Int64,Int64,Int64}[(1,1,1)])  )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = JAC.give("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-PhotoExcitationFluores-new.sum")
    success = true
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoExcitationFluores-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoExcitationFluores-new.sum"),  "xxx", 10) 
    testPrint("testModule_PhotoExcitationFluores()::", success)
    return(success)  
end



"""
`JAC.testModule_PhotoIonization(; short::Bool=true)`  ... tests on module JAC.PhotoIonization.
"""
function testModule_PhotoIonization(; short::Bool=true) 
    JAC.define("print summary: open", "test-PhotoIonization-new.sum")
    printstyled("\n\nTest the module  PhotoIonization  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(36.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 800),
                            initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                            process = JAC.Photo, 
                            processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 4000.], false, true, true, true, 
                                            true, Tuple{Int64,Int64}[(1,1), (1,2)], ExpStokes(1., 0., 0.)) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoIonization-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoIonization-new.sum"), "Total photoionization c", 8) 
    testPrint("testModule_PhotoIonization()::", success)
    return(success)  
end



"""
`JAC.testModule_PhotoRecombination(; short::Bool=true)`  ... tests on module JAC.PhotoRecombination.
"""
function testModule_PhotoRecombination(; short::Bool=true) 
    JAC.define("print summary: open", "test-PhotoRecombination-new.sum")
    printstyled("\n\nTest the module  PhotoRecombination  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(36.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600),
                            initialConfigs=[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^6")], 
                            process = JAC.Rec, 
                            processSettings=PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 20.], [0.], 
                                            false, true, true, true, true, Tuple{Int64,Int64}[(1,1)]) )
    wb = perform(wa)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-PhotoRecombination-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-PhotoRecombination-new.sum"), "Photorecombination cross sections for ", 12) 
    testPrint("testModule_PhotoRecombination()::", success)
    return(success)  
end



"""
`JAC.testModule_Radiative(; short::Bool=true)`  ... tests on module JAC.Radiative.
"""
function testModule_Radiative(; short::Bool=true) 
    JAC.define("print summary: open", "test-Radiative-new.sum")
    printstyled("\n\nTest the module  Radiative  ... \n", color=:cyan)
    ### Make the tests
    wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                            initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")], 
                            process = JAC.RadiativeX, 
                            processSettings=Radiative.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, 
                            true, Tuple{Int64,Int64}[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)], 0., 0., 10000. ) )
    ##x streamDummy = open(pwd() * "/runtests.dummy", "w")
    ##x redirect_stdout(streamDummy) do   
    wb = perform(wa)   
    ##x end
    ##x close(streamDummy)
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-Radiative-approved.sum"), 
                                joinpath(@__DIR__, "..", "test", "test-Radiative-new.sum"), "Einstein coefficients, t", 100) 
    testPrint("testModule_Radiative()::", success)
    return(success)  
end



"""
`JAC.testModule_RayleighCompton(; short::Bool=true)`  ... tests on module JAC.RayleighCompton.
"""
function testModule_RayleighCompton(; short::Bool=true) 
    JAC.define("print summary: open", "test-RayleighCompton-new.sum")
    printstyled("\n\nTest the module  RayleighCompton  ... \n", color=:cyan)
    ### Make the tests
    ###
    JAC.define("print summary: close", "")
    # Make the comparison with approved data
    printTest, iostream = JAC.give("test flag/stream")
    println(iostream, "Make the comparison with approved data for ... test-RayleighCompton-new.sum")
    success = true
    ## success = testCompareFiles( joinpath(@__DIR__, "..", "test", "approved", "test-RayleighCompton-approved.sum"), 
    ##                             joinpath(@__DIR__, "..", "test", "test-RayleighCompton-new.sum"), "xxx", 100) 
    testPrint("testModule_RayleighCompton()::", success)
    return(success)  
end


