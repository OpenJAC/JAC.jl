
using    GSL, HypergeometricFunctions, Plots, Printf, SpecialFunctions, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..Radial, ..ManyElectron, 
            ..Nuclear, ..Pulse, ..TableStrings

"""
`StrongField.computeSphericalAmplitudesUncoupled(comp::StrongField.Computation)`  
    ... computes all necessary spherical amplitudes for the given initial and final levels, the Volkov states
        and polarization and by using the pulse envelope integrals. A newAmplitudes::Array{SphericalAmplitude,1} is returned.
"""
function  computeSphericalAmplitudesUncoupled(comp::StrongField.Computation)
    
    #Extract the initial orbital of the active electron from the many-electron comp.initialLevel and set quantum numbers
    initialOrbitals = comp.initialLevel.basis.orbitals

    #Find highest lying orbital (smallest ionization potential)
    defaultSubshell = [sh for (sh,or) in initialOrbitals][1] #This is not nice; must be a better way to simply get a default element from a Dict
    initialOrbital = initialOrbitals[defaultSubshell]
    minIonPotential = abs(initialOrbital.energy)
    for (subshell,orbital) in initialOrbitals
        if   abs(orbital.energy) < minIonPotential
            initialOrbital = orbital
            minIonPotential = abs(orbital.energy)
        end
    end
    
    initialEneV = convertUnits("energy: from atomic to eV", initialOrbital.energy)
    println("")
    println("Selected orbital $(initialOrbital.subshell) as active electron for the strong-field ionization process. Binding energy: $(initialEneV) eV.")
    println("")
    
    ls = LevelSymmetry(initialOrbital.subshell)
    n = initialOrbital.subshell.n;      l = (ls.J.num+1)/2;    initialEn = initialOrbital.energy
    if  (sign((-1)^l) == -1 && ls.parity == plus::Parity) || (sign((-1)^l) == 1 && ls.parity == minus::Parity)
        l = l - 1
    end
    l = floor(Int,l)
    
    if comp.settings.hydrogenic1s
        n = 1
        l = 0
    end
    
    lpMin = abs(l-1)
    lpMax = l+1
    
    # First determine which spherical SFA amplitudes need to be computed before the actual computation starts
    mNum = 1    #If comp.Average == false
    miAverageFactor = 1.      #If comp.Average == false
    if( comp.settings.mAverage )
        mNum = 2*l + 1
        miAverageFactor = 1.0/sqrt( mNum )
    end
    sfaAmplitudes = StrongField.determineSphericalAmplitudes(comp.observable,mNum,1)
    
    # Compute the requested SFA amplitudes
    for jm = 1:mNum
        m = -l + (jm-1) 
        for  jAmp = 1:length(sfaAmplitudes[jm,1,:]) #amp in  sfaAmplitudes[jm,1]
            amp = sfaAmplitudes[jm,1,jAmp]
            thetap        = amp.theta;   phip = amp.phi;     energyp = amp.energy
            envIntegrals  = StrongField.computeEnvelopeIntegrals(comp.envelope, comp.beam, comp.polarization, thetap, phip, energyp, initialEn)
            fVolkovPlus, fVolkovMinus, fVolkovSquared = envIntegrals
            #
            reducedMEArray = ComplexF64[]
            for lp = lpMin:lpMax
                push!( reducedMEArray, pReducedMEHydrogenicUncoupled(energyp, lp, n, l, initialEn, comp.volkov) )
            end
                
            wminus = 0.0im;    wplus = 0.0im
            # Collect contributions from all l_p, q terms
            for  lp = lpMin:lpMax
                reducedME = reducedMEArray[lp-lpMin+1]
                
                KWignerEckart = 1.0/sqrt(2*lp+1) #Factor in the Wigner-Eckart theorem (CONVENTION that needs to fit the reduced ME!)
                for  q in [-1, 0, 1]
                    wplus = wplus + KWignerEckart * AngularMomentum.sphericalYlm(lp, m-q, amp.theta, amp.phi) * (-1)^q  * 
                            Basics.determinePolarizationVector(q, comp.polarization, star=false)         *
                            AngularMomentum.ClebschGordan_old( AngularJ64(Rational(l)), AngularM64(Rational(m)), AngularJ64(1), AngularM64(Rational(-q)), AngularJ64(Rational(lp)), AngularM64(Rational(m-q)) ) * reducedME
                            
                    wminus = wminus + KWignerEckart * AngularMomentum.sphericalYlm(lp, m+q, amp.theta, amp.phi)           * 
                            Basics.determinePolarizationVector(q, comp.polarization, star=true)          *
                            AngularMomentum.ClebschGordan_old( AngularJ64(Rational(l)), AngularM64(Rational(m)), AngularJ64(1), AngularM64(Rational(q)), AngularJ64(Rational(lp)), AngularM64(Rational(m+q)) ) * reducedME
                end
            end
            
            scalarProd = scalarProdBoundContHydrogenicUncoupled(energyp, l, n, l, m, initialEn, comp.volkov)
        
            wa = (fVolkovPlus * wminus + fVolkovMinus * wplus) + fVolkovSquared * AngularMomentum.sphericalYlm(l, m, amp.theta, amp.phi) * scalarProd
            wa = - im/sqrt(2*pi) * wa
            wa = miAverageFactor * wa
            
            sfaAmplitudes[jm,1,jAmp] = SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa)
        
            if( comp.settings.printAmplitudes )
                println(">> $(SphericalAmplitude(amp.energy, amp.theta, amp.phi, wa))")
            end
        end #for amp
    end #for jmj
            
    return( sfaAmplitudes )
end

