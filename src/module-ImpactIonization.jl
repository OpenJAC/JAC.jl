
"""
`module  JAC.ImpactIonization`  
    ... a submodel of JAC that contains all methods for computing electron impact ionization cross sections.
"""
module ImpactIonization

    using  Printf, ..Basics, ..Cubature, ..ManyElectron, ..Nuclear, ..TableStrings, ..Defaults, ..Radial, ..RadialIntegrals
    
    """
    `abstract type ImpactIonization.AbstractApproximation` 
        ... defines an abstract and a number of singleton types for the computation of empirical EII cross sections in different 
            approximations.
        
        + struct BEBapproximation        
            ... to apply a modified BEB approximation for non-relativistic impact energies.
        + struct BEDapproximation        
            ... to apply a modified BED approximation for non-relativistic impact energies.  
        + struct RelativisticBEBapproximation        
            ... to apply a modified relativistic BEB approximation for non-relativistic & relativistic impact energies.
        + struct RelativisticBEDapproximation
            ... to apply a modified relativistic BEB approximation due to Uddin et al. (20xx) for non-relativistic & 
                relativistic impact energies.
        + struct FittedBEDapproximation  #  MUIBEDapproximation   
            ... to apply a modified relativistic BEB approximation due to Uddin et al. (20xx) for non-relativistic & 
                relativistic impact energies, and by using fitting coefficients for K, L, and M subshells. 
    """
    abstract type  AbstractApproximation  
    end
    
    struct         BEBapproximation               <:  AbstractApproximation   end
    struct         BEDapproximation               <:  AbstractApproximation   end
    struct         RelativisticBEBapproximation   <:  AbstractApproximation   end
    struct         RelativisticBEDapproximation   <:  AbstractApproximation   end
    struct         FittedBEDapproximation         <:  AbstractApproximation   end
    
    export BEBapproximation, BEDapproximation, RelativisticBEBapproximation, RelativisticBEDapproximation, FittedBEDapproximation
    

    """
    `struct  ImpactIonization.Settings  <:  AbstractEmpiricalSettings`  
        ... defines a type for the settings of empirical electron-impact ionization cross sections.
        
        + approximation              ::ImpactIonization.AbstractApproximation 
            ... approximation for computing empirical EII cross sections.
        + impactEnergies             ::Array{Float64,1}                       
            ... List of impact-energies of the incoming elecgtrons.
        + calcDifferentialCs         ::Bool                                   
            ... True if energy-differential cross sections are to be calculated, and false otherwise.
        + calcPartialCs              ::Bool                                   
            ... True if partial (shell-dependent) cross sections are to be calculated, and false otherwise.
        + calcTotalCs                ::Bool                                   
            ... True if total cross sections are to be calculated, and false otherwise.
        + shellSelection             ::ShellSelection      ... Describe the selected shells, if any.                               ... 
    """    
    struct Settings   <:  AbstractEmpiricalSettings
        approximation                ::ImpactIonization.AbstractApproximation
        impactEnergies               ::Array{Float64,1}
        calcDifferentialCs           ::Bool 
        calcPartialCs                ::Bool 
        calcTotalCs                  ::Bool 
        shellSelection               ::ShellSelection
    end      
    

    """
    `JAC.ImpactIonization.Settings()`  ... constructor for the default values of empirical electron-impact ionization cross sections.
    """
    function Settings()
       Settings(ImpactIonization.AbstractApproximation(), Float64[], Float64[], false, false, false, ShellSelection(), Shell())
    end
    
    
    ## `Base.show(io::IO, settings::ImpactIonization.Settings)`  
    ##  ... prepares a proper printout of the variable settings::ImpactIonization.Settings.
    function Base.show(io::IO, settings::ImpactIonization.Settings) 
        println(io, "approximation:           $(settings.approximation)  ")
        println(io, "impactEnergies:          $(settings.impactEnergies)  ")
        println(io, "calcPartialCs:           $(settings.calcPartialCs)  ")
        println(io, "calcTotalCs:             $(settings.calcTotalCs)  ")
        println(io, "shellSelection:          $(settings.shellSelection)  ")     
    end


    """
    `struct  ImpactIonization.CrossSection`  
        ... defines a type for dealing with (shell-dependent) partial EII cross sections that are based on empirical computations.

        + subshell          ::Subshell          ... subshell for which the (differential and partial) cross sections are calculated.
        + impactEnergy      ::Float64           ... energy of impact electron [Hartree].
        + partialCS         ::Float64           ... partial cross sections for this shell.
        + electronEnergies  ::Array{Float64,1}  ... energies of the ionized electron.
        + differentialCS    ::Array{Float64,1}  ... differential cross sections for this shell.
    """
    struct  CrossSection
        subshell            ::Subshell
        impactEnergy        ::Float64
        partialCS           ::Float64
        electronEnergies    ::Array{Float64,1}
        differentialCS      ::Array{Float64,1}
    end 


    # `Base.show(io::IO, cs::ImpactIonization.CrossSection)`  ... prepares a proper printout of the variable line::ImpactIonization.CrossSection.
    function Base.show(io::IO, cs::ImpactIonization.CrossSection) 
        println(io, "subshell:           $(cs.subshell)  ")
        println(io, "impactEnergy:       $(cs.impactEnergy)  ")
        println(io, "partialCS:          $(cs.partialCS)  ")
        println(io, "electronEnergies:   $(cs.electronEnergies)  ")
        println(io, "differentialCS:     $(cs.differentialCS)  ")
    end

   
        
    """
    `ImpactIonization.computeCrossSections(basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings; output::Bool=true)`  
        ... to compute the EII cross sections for an atom or ion in BEBapproximation. While the basis describes the
            configuration of the atom, the settings provide access to the selected subshells and energies; 
            a list of crossSections::Array{ImpactIonization.CrossSection,1} is returned.
    """    
    function  computeCrossSections(basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings; output::Bool=true)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        println("")
        printstyled("ImpactIonization.computeCrossSections(): Use approximation $(settings.approximation) ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        newCrossSections = ImpactIonization.CrossSection[]
        crossSections    = ImpactIonization.determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)
        # 
        for  cs in crossSections
            newCs = ImpactIonization.computeCrossSections(settings.approximation, cs, basis, nm, settings) 
            push!( newCrossSections, newCs)
        end
        # Print all results to screen
        ImpactIonization.displayCrossSections(stdout, newCrossSections, settings)
        if  printSummary   ImpactIonization.displayCrossSections(iostream, newCrossSections, settings)     end
        #
        if    output    return( newCrossSections )
        else            return( nothing )
        end
    end
   
        
    """
    `ImpactIonization.computeCrossSections(approx::BEBapproximation, cs::ImpactIonization.CrossSection,
                                           basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)`  
        ... to compute the particular EII cross sections (cs) for an atom or ion in BEBapproximation. While the basis describes the
            configuration of the atom, the settings provide access to the selected subshells and energies; 
            a cs::ImpactIonization.CrossSection is returned.
    """    
    function  computeCrossSections(approx::BEBapproximation, cs::ImpactIonization.CrossSection,
                                   basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)
        
        #== Yuan-Cheng: only use those lines which you need in the BEBapproximation
        # Determine the binding energy, kinetic energy and occupation number for the subshell
        bindingEnergy = - basis.orbitals[subshell].energy
        a = basis.orbitals[subshell]
        kineticEnergy = RadialIntegrals.isotope_nms(a, a, 0.0, grid)
        occupationNumber = Basics.computeMeanSubshellOccupation(subshell, basis)
        eC = ImpactIonization.effectiveChargeNumber(nm, subshell, subshellList, basis)
        # Transform the imput energies into atomic unit and get the default values
        iE = Defaults.convertUnits("energy: from eV to atomic", impactEnergy)
        mC2 = Defaults.getDefaults("mc^2") 
        alpha = Defaults.getDefaults("alpha")
        R = convertUnits("energy: from atomic to eV", 0.5)
        # Define parameters using symbols in BEB & BED formulas
        t = iE / bindingEnergy;     u = kineticEnergy / bindingEnergy;  
        q = occupationNumber 
        # Define parameters using symbols in relativisticBEB & BED formulas
        tP = iE / mC2;  uP = kineticEnergy / mC2;  bP = bindingEnergy / mC2
        relTpra = 1 - 1/(1 + tP)^2; relBpra = 1 - 1/(1 + bP)^2; relUpra = 1 - 1/(1 + uP)^2
        # Calculate partial ionization cross section by approximations.
        partialCS = 0.0 
        ==#
        bindingEnergy = q = t = u = 0.0
        partialCS     = 0.0 
        eC            = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
        @show cs.subshell, eC
        #
        # Calculate the partial and total cross sections
        csFactor = pi * q * (1/bindingEnergy)^2/(t + (u + 1)/(eC+1))
        wc       = 1 - 1/t - log(t)/(1+t)
        wd       = log(t) * 1/2 * (1 - 1/t^2 )            
        if t < 1   else    partialCS = csFactor * (wc + wd)     end 
            
        newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
        
        return( newCs )
    end
   
        
    """
    `ImpactIonization.computeCrossSections(approx::BEDapproximation, cs::ImpactIonization.CrossSection,
                                           basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)`  
        ... to compute the EII cross sections (cs) for an atom or ion in BEDapproximation. While the basis describes the
            configuration of the atom, the settings provide access to the selected subshells and energies; 
            a cs::ImpactIonization.CrossSection is returned.
    """    
    function  computeCrossSections(approx::BEDapproximation, cs::ImpactIonization.CrossSection,
                                   basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)
        
        # Yuan-Cheng: only use those lines which you need in the BEDapproximation
           
        newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, Float64[], Float64[], partialCS)
        
        return( newCs )
    end
   
        
    """
    `ImpactIonization.computeCrossSections(approx::FittedBEDapproximation, cs::ImpactIonization.CrossSection,
                                           basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)`  
        ... to compute the EII cross sections (cs) for an atom or ion in FittedBEDapproximation. While the basis describes the
            configuration of the atom, the settings provide access to the selected subshells and energies; 
            a cs::ImpactIonization.CrossSection is returned.
    """    
    function  computeCrossSections(approx::FittedBEDapproximation, cs::ImpactIonization.CrossSection,
                                   basis::Basis, nm::Nuclear.Model, settings::ImpactIonization.Settings)
        #
        function scalingFactor(subshell::Subshell, nm::Nuclear.Model)
            # Determine a scaling factor in the FittedBEDapproximation
            parameters = [ (Subshell("1s_1/2"),  0.,    0.), 
                           (Subshell("2s_1/2"), -0.02,  0.001),  (Subshell("2p_1/2"), 0.13, 0.11),  (Subshell("2p_3/2"), 0., 0.), 
                           (Subshell("2s_1/2"), 0., 0.),  (Subshell("3p_1/2"), 0., 0.), (Subshell("3p_3/2"), 0., 0.), 
                           (Subshell("3d_3/2"), 0., 0.),  (Subshell("3d_5/2"), 0., 0.)  ]
            #== if   subshell.n == 1                          m = 0.0;   lambda = 0.0
            elseif  subshell.n == 2      
                if      subshell.kappa == -1              m =-0.02;  lambda = 0.001
                elseif  subshell.kappa ==  1              m = 0.13;  lambda = 0.11
                elseif  subshell.kappa == -2              m = 0.04;  lambda = 0.32
                end
            elseif  subshell.n == 3      
                if      subshell.kappa == -1              m =-0.01;  lambda = 0.01
                elseif  subshell.kappa ==  1              m =-0.07;  lambda = 0.01
                elseif  subshell.kappa == -2              m = 0.05;  lambda = 0.04
                elseif  subshell.kappa ==  2              m = 0.45;  lambda = 0.05
                elseif  subshell.kappa == -3              m = 0.40;  lambda = 0.01
                end
            else    println("Can only calculate K- L- and M-shells!");   error("stop a")
            end ==#
            hasfound = false
            for  par in parameters   
                if par[1] == subshell  m = par[2];  lambda = par[3];   hasfound = true;    break    end
            end
            if  !hasfound   error("Can only calculate K- L- and M-shells!")     end
                
            scalingFactor = 1 + m * nm.Z^lambda
            return( scalingFactor )
        end
        
        # Yuan-Cheng: only use those lines which you need in the FittedBEDapproximation
           
        newCs = ImpactIonization.CrossSection(cs.subshell, Float64[], Float64[], partialCS)
        
        return( newCs )
    end
    
    
    # Yuan-Cheng: ... add the procedures for the other approximations here !!
    
    
     
    """
    `ImpactIonization.computeCrossSections(subshell::Subshell, impactEnergy::Float64, basis::Basis, 
                                           nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings, 
                                           subshellList::Array{Subshell,1})`  
        ... OLD: to compute the EII cross sections for given subshell and impact energy; a crossSection::Float64 is returned.
            OLD: to be deleted.
    """    
    function  computeCrossSections_old(subshell::Subshell, impactEnergy::Float64, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings, subshellList::Array{Subshell,1})
        # Determine the binding energy, kinetic energy and occupation number for the subshell
        bindingEnergy = - basis.orbitals[subshell].energy
        a = basis.orbitals[subshell]
        kineticEnergy = RadialIntegrals.isotope_nms(a, a, 0.0, grid)
        occupationNumber = Basics.computeMeanSubshellOccupation(subshell, basis)
        eC = ImpactIonization.effectiveChargeNumber(nm, subshell, subshellList, basis)
        # Transform the imput energies into atomic unit and get the default values
        iE = Defaults.convertUnits("energy: from eV to atomic", impactEnergy)
        mC2 = Defaults.getDefaults("mc^2") 
        alpha = Defaults.getDefaults("alpha")
        R = convertUnits("energy: from atomic to eV", 0.5)
        # Define parameters using symbols in BEB & BED formulas
        t = iE / bindingEnergy;  u = kineticEnergy / bindingEnergy;  
        q = occupationNumber 
        # Define parameters using symbols in relativisticBEB & BED formulas
        tP = iE / mC2;  uP = kineticEnergy / mC2;  bP = bindingEnergy / mC2
        relTpra = 1 - 1/(1 + tP)^2; relBpra = 1 - 1/(1 + bP)^2; relUpra = 1 - 1/(1 + uP)^2
        # Calculate partial ionization cross section by approximations.
        partialCS = 0.0         
        if settings.approximation == ImpactIonization.RelativisticBEBapproximation()
            csFactor = 4 * pi * alpha^4 * q /( 2 * bP) * 1 /(relTpra + (relBpra + relUpra)/(eC +1))
            wa = log( relTpra / (1 - relTpra) ) - relTpra - log(2 * bP) 
            wb = 1 - 1/t - log(t)/(1 + t) * (1 + 2 * tP)/(1 + tP / 2)^2 + bP^2/(1 + tP / 2)^2 *(t - 1)/2        
            analyticDipole = 0.5 * (1 - 1 / t^2) 
            if t < 1   partialCS = 0.0
            else   partialCS = csFactor * (wa * analyticDipole +  wb)
            end
        elseif settings.approximation == ImpactIonization.RelativisticBEDapproximation()
            Hfactor = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
            SRfactor = 4 * pi * alpha^2 * q /(relTpra + (relBpra + relUpra)/(eC +1))
            FRfactor = 64/alpha * q * (sqrt(relBpra))^3 /relTpra
            k = 2*bindingEnergy
            f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
            xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+tP), iE-bindingEnergy) 
            (val,err)=hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
            Gfactor = val
            if t < 1   partialCS = 0.0
            else partialCS = SRfactor *Hfactor + FRfactor *Gfactor
            end
        elseif settings.approximation == BEBapproximation()
            csFactor = pi * q * (1/bindingEnergy)^2/(t + (u + 1)/(eC+1))
            wc = 1 - 1/t - log(t)/(1+t)
            wd = log(t) * 1/2 * (1 - 1/t^2 )            
            if t < 1   partialCS = 0.0
            else partialCS = csFactor * (wc + wd)
            end  
        elseif settings.approximation == BEDapproximation()
            Hfactor = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
            Sfactor = 4 * pi * q /(iE +(bindingEnergy + kineticEnergy)/(eC + 1))
            Ffactor = 32 * q * (sqrt(2*bindingEnergy))^3 /iE 
            k = 2*bindingEnergy
            f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
            xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+tP), iE-bindingEnergy) 
            (val,err)=hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
            Gfactor = val
            if t < 1   partialCS = 0.0
            else partialCS = Sfactor *Hfactor + Ffactor *Gfactor
            end
        elseif settings.approximation == ImpactIonization.FittedBEDapproximation()
            scalingFactor = ImpactIonization.scalingFactor(subshell, nm)
            Hfactor = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
            Sfactor = 4 * pi * q /(iE +(bindingEnergy + kineticEnergy)/(eC + 1))
            Ffactor = 32 * q * (sqrt(2*bindingEnergy))^3 /iE 
            k=2*bindingEnergy
            f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
            xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+tP), iE-bindingEnergy)  
            (val,err)=hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
            Gfactor = val
            RfactorB = (((1 + t)*(t + 2/bP)*(1 + 1/bP)^2)/((1 + 2/bP)/(bP^2) + t*(t + 2/bP)*(1 + 1/bP)^2))^1.5
            RfactorLM = (1 + 2/bP)/(0.6 * t +2/bP) * ((t + 1/bP)/(1 + 1/bP))^2.11
            partialCS = scalingFactor * (RfactorLM *Sfactor *Hfactor + RfactorB *Ffactor *Gfactor)
        end 
        return( partialCS )
    end


    """
    `ImpactIonization.determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)`  
        ... to determine an array of EII cross sections for all requested shells of the given multiplett that need to be calculated in these
            computations; an Array{ImpactIonization.CrossSection,1} is returned.
    """
    function  determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)
        crossSections = ImpactIonization.CrossSection[]
        if  settings.shellSelection.active
            subshells = Basics.generateSubshellList(settings.shellSelection.shells) 
        else
            confs     = Basics.extractNonrelativisticConfigurations(basis)
            shells    = Basics.extractShellList(confs)
            subshells = Basics.generateSubshellList(shells)
        end
        # Generate the list of cross sections
        for  subsh in subshells
            for  iE in  settings.impactEnergies
                push!(crossSections, ImpactIonization.CrossSection(subsh, iE, 0., Float64[], Float64[]))
            end
        end
        
        return( crossSections )
    end


    """
    `ImpactIonization.displayCrossSections(stream::IO, crossSections::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)`  
        ... displays the EII cross sections in a neat tabular form; nothing is returned.
    """
    function  displayCrossSections(stream::IO, crossSections::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)
        # Display the partial (shell-dependent) cross sections
        if  settings.calcPartialCs 
            nx = 120 
            println(" ")
            println("  Partial ionization cross sections for atoms by electron impact:")
            println(" ")
            println("  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(18, "subshell"         ; na=0);                       
            sb = sb * TableStrings.hBlank(18)
            sa = sa * TableStrings.center(18, "binding energy"   ; na=2);                       
            sb = sb * TableStrings.center(18, TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(18, "kinetic energy"   ; na=4)               
            sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=0)
            sa = sa * TableStrings.center(12, "occupation number"  ; na=4)             
            sb = sb * TableStrings.hBlank(30) 
            sa = sa * TableStrings.center(12, "impact energies"  ; na=4)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)            
            sa = sa * TableStrings.center(28, "Cross section"; na=3)      
            sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section"); na=3)
            println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
            #
            for  cs  in crossSections
                @show cs.subshell, cs.impactEnergy, cs.partialCS
                ## Yuan-Cheng ... Now extract all data from cs only !!
                ##
            end
            println("  ", TableStrings.hLine(nx))    
        end
        #
        # Display the total EII cross sections
        if  settings.calcTotalCs 
            nx = 120 
            println(" ")
            println("  Total ionization cross sections for atoms by electron impact:")
            println(" ")
            println("  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            ## Yuan-Cheng ... Now display the total cross sections, either from the cs data ... or by further 
            ##                summation over these data
            println("  ", TableStrings.hLine(nx))    
        end
        #
        return( nothing )
    end
    

    """
    `ImpactIonization.displayCrossSections(basis::Basis, grid::Radial.Grid, nm::Nuclear.Model, impactEnergies::Array{Float64,1}, shells::Array{Shell,1}, settings::ImpactIonization.Settings, configs::Array{Configuration,1})`  
        OLD: ... displays the EII cross sections in a neat tabular form; nothing is returned.
        OLD: ... to be deleted.
    """
    function  displayCrossSections_old(basis::Basis, grid::Radial.Grid, nm::Nuclear.Model, impactEnergies::Array{Float64,1}, shells::Array{Shell,1}, settings::ImpactIonization.Settings, configs::Array{Configuration,1})
        subshellList = Basics.extractRelativisticSubshellList(configs)
        if  settings.calcPartialCs 
            nx = 120 
            println(" ")
            println("  Partial ionization cross sections for atoms by electron impact:")
            println(" ")
            println("  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(18, "subshell"         ; na=0);                       
            sb = sb * TableStrings.hBlank(18)
            sa = sa * TableStrings.center(18, "binding energy"   ; na=2);                       
            sb = sb * TableStrings.center(18, TableStrings.inUnits("energy"); na=4)
            sa = sa * TableStrings.center(18, "kinetic energy"   ; na=4)               
            sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=0)
            sa = sa * TableStrings.center(12, "occupation number"  ; na=4)             
            sb = sb * TableStrings.hBlank(30) 
            sa = sa * TableStrings.center(12, "impact energies"  ; na=4)             
            sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)            
            sa = sa * TableStrings.center(28, "Cross section"; na=3)      
            sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section"); na=3)
            println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
            subshells = ImpactIonization.determineCrossSections(basis, shells, settings)
            for subshell in subshells
                a = basis.orbitals[subshell]
                bindingEnergy = - basis.orbitals[subshell].energy
                kineticEnergy = 0.5* RadialIntegrals.isotope_nms(a, a, 0.0, grid)
                occupationNumber = Basics.computeMeanSubshellOccupation(subshell, basis)
                sa  = ""; 
                sa = sa * TableStrings.center(18, string(subshell); na=6)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to eV", bindingEnergy)) * "        "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to eV", kineticEnergy)) * "           "
                sa = sa * @sprintf("%.6e", occupationNumber) * "   "
                println(sa)
                for impactEnergy in impactEnergies
                    sb   = TableStrings.hBlank(86)
                    sb = sb * @sprintf("%.6e", impactEnergy) * "             "    
                    wa = ImpactIonization.computeCrossSections(subshell, impactEnergy, basis, nm, grid, settings, subshellList)
                    sb = sb * @sprintf("%.10e", Defaults.convertUnits("cross section: from atomic to barn", wa))
                    println(sb)
                end 
            end
        end
        
        println("  ", TableStrings.hLine(nx))    
        return( nothing )
    end


    """
    `ImpactIonization.effectiveCharge(subshell::Subshell, nm::Nuclear.Model, basis::Basis)`  
        ... determined the effective charge that is felt by an electron in the given subshell;
            a zEff::Float64 is returned.
    """
    function  effectiveCharge(subshell::Subshell, nm::Nuclear.Model, basis::Basis)
        wz        = nm.Z
        confs     = Basics.extractNonrelativisticConfigurations(basis)
        shells    = Basics.extractShellList(confs)
        subshells = Basics.generateSubshellList(shells)
        enSub     = basis.orbitals[subshell].energy
        for  subsh  in subshells   
            if  basis.orbitals[subsh].energy < enSub
                wz = wz - Basics.computeMeanSubshellOccupation(subsh, basis) 
            end 
        end
        
        return(wz)
    end

end # module
