using DelimitedFiles, Plots


"""
`StrongField.plot( comps::Array{StrongField.Computation,1}, results::Array{Dict{String,Any},1}, energyScale::String = "atomic", probabilityScaling::String = "linear", dataLabel::String = "StrongFieldData" )`  
    ... generates a graphical representation of the observable (StrongField.SfaEnergyDistribution, StrongField.SfaMomentumDistribution, StrongField.SfaAzimuthalAngularDistribution or StrongField.SfaPolarAngularDistribution) with the results of all computations comps combined in one plot. 
    - All comps.observable need to be equal!
    - The probabilities for all comps need to be at the same grid points (either energies, momenta or angles)
    - energyScale determines the scaling of the energy axis either in atomic units (energyScale = "atomic") or in units of hbar*omega (energyScale = "omega").
    - The y-axis is scaled either linearly (probabilityScaling = "linear") or logarithmically (probabilityScaling = "log")
"""
function plot( comps::Array{StrongField.Computation,1}, results::Array{Dict{String,Any},1}, energyScale::String = "atomic", probabilityScaling::String = "linear", dataLabel::String = "StrongFieldData" )

    #---Check if all observables are equal---
    observable = typeof(comps[1].observable)
    for comp in comps
        if  typeof(comp.observable) != observable
            error("Unequal types of observables for data given to StrongField.plot().")
            return nothing
        end
    end

    #Define colors and linestyles for plot of SfaEnergyDistribution and SfaAngularDistribution
    colors = [:black,:blue,:red,:green,:purple]
    styles = filter((s->begin
            s in Plots.supported_styles()
        end), [:solid, :dash, :dot, :dashdot, :dashdotdot])

    #---Photoelectron energy spectra---
    if  observable == StrongField.SfaEnergyDistribution  
    
        for k = 1:size(comps,1) #Loop over all comps
        
            #prepare data and rescale the x-axis if neccessary
            if  energyScale == "atomic"
                energyLabel = "E (a.u.)"
                energies = results[k]["energy distribution"].energies
            elseif  energyScale == "omega"
                energyLabel = "E/w"
                energies = results[k]["energy distribution"].energies / comps[k].beam.omega
            end
            probabilities = results[k]["energy distribution"].probabilities
            gr() #sets the plotting backend to the package "GR"
            
            
            #set scaling of the y-axis
            scaling = :identity
            if  probabilityScaling == "log"
                        scaling = :log10
            end
            
            #generate the plot
            if  k == 1 #initial plot generation
                Plots.plot(energies, probabilities, 
                            title = "Photoelectron energy spectrum",
                            xlabel = energyLabel,
                            ylabel = "P(E)",
                            xscale = :identity,
                            yscale = scaling,
                            framestyle = :box,
                            legend = :none,
                            #markershape = :circle,
                            #markershape = :none,
                            line = (2,styles[1]),
                            linecolor = colors[1],
                            gridlinewidth = 2,
                            tickfontsize = 10,
                            labelfontsize = 10,
                            labelfontfamily = "Latin Modern Roman",
                            titlefontfamily = "Latin Modern Roman"
                            )
            else    #add to the plot already generated for k = 1
                Plots.plot!(energies, probabilities, 
                            title = "Photoelectron energy spectrum",
                            xlabel = energyLabel,
                            ylabel = "P(E)",
                            xscale = :identity,
                            yscale = scaling,
                            framestyle = :box,
                            legend = :none,
                            #markershape = :circle,
                            #markershape = :none,
                            line = (2,styles[k]),
                            linecolor = colors[k],
                            gridlinewidth = 2,
                            tickfontsize = 10,
                            labelfontsize = 10,
                            labelfontfamily = "Latin Modern Roman",
                            titlefontfamily = "Latin Modern Roman"
                            )
            end
                            
        end
                            
        #export the plot as png-file
        savefig(dataLabel * "-energy_spectrum.pdf")
    #-----------------------------------
    
    #---Photoelectron momentum distributions - DOES NOT WORK YET---
    elseif  observable == StrongField.SfaMomentumDistribution
            #prepare data
            #pList = results["momentum distribution"].momenta
            #probList = results["momentum distribution"].probabilities

            #generate the plot
            #pyplot()
            #r = range(0,stop=10,length=11)
            #theta = range(0,stop=360,length=361)
            #f(r,theta) = r^2
            #println("$(f.(r,theta'))")
            #Plots.plot( heatmap( f.(r,theta'), proj=:polar ) )
            #println("$(transpose(probList))")
            #Plots.plot( heatmap( pList, transpose(probList), proj=:polar ) )
            #imshow(pList,transpose(probList))
            #probList = probList[2:-1]
            #matplotlib.pyplot.pcolormesh(pList[:,1], pList[:,2], transpose(probList))
            contourf( pList[:,1], pList[:,2], probList )
            #matplotlib.pyplot.imshow(pList)
            #heatmap( pList[:,1], pList[:,2], probList, proj=:polar, legend=true
                            #heatmap( probList
            #                )

            #p = Plots.plot( #heatmap( pList[:,1], pList[:,2], probList
            #               Plots.GR.polarheatmap( probList 
                        ##heatmap( probList
            #                 )
            #              )
            
            
            #export the plot as png-file
            #savefig(dataLabel * "-momentum_distribution.pdf")
    #-----------------------------------------
    
    #---Photoelectron angular distributions---
    elseif  observable == StrongField.SfaAzimuthalAngularDistribution || observable == StrongField.SfaPolarAngularDistribution
        gr() #sets the plotting backend to the package "GR"
            
        for k = 1:size(comps,1) #Loop over all comps
            
            #prepare data
            if  observable == StrongField.SfaAzimuthalAngularDistribution
                angleList = results[k]["angular distribution"].phis
                plotTitle = "Photoelectron angular distribution (azimuthal)"
            elseif observable == StrongField.SfaPolarAngularDistribution
                angleList = results[k]["angular distribution"].thetas
                plotTitle = "Photoelectron angular distribution (polar)"
            end
            probList = results[k]["angular distribution"].probabilities
            
            
            lineLabel = ""
            if energyScale == "atomic"
                lineLabel = string( round( results[k]["angular distribution"].energy, digits = 2) ) * " a.u."
                legendTitle = "Energy"
            elseif energyScale == "omega"
                lineLabel = string( round( results[k]["angular distribution"].energy / comps[k].beam.omega, digits = 1 ) ) * " w"
                legendTitle = "Energy"
            end
            
            if  k == 1 #initial plot generation
                GR.polar(angleList, probList)
                Plots.plot(angleList,probList,
                            title = plotTitle,
                            proj = :polar,
                            line = (2,styles[1]),
                            linecolor = colors[1], 
                            labelfontfamily = "Latin Modern Roman",
                            titlefontfamily = "Latin Modern Roman",
                            gridlinewidth = 2,
                            tickfontsize = 10,
                            labelfontsize = 10,
                            label = lineLabel, 
                            legendtitle = legendTitle,
                            legend = true
                            
                            
                            )
            else
                Plots.plot!(angleList,probList,
                            proj = :polar,
                            line = (2,styles[k]),
                            linecolor = colors[k],
                            gridlinewidth = 2,
                            tickfontsize = 10,
                            labelfontsize = 10,
                            label = lineLabel,
                            labelfontfamily = "Latin Modern Roman",
                            
                            )
            end
        
        end
        
        #export the plot as pdf-file
        savefig(dataLabel * "-angular_distribution.pdf")
    #-----------------------------------------
            
    #---Not a valid obserable---
    else     
            error("Undefined observable for strong-field computations in StrongField.plot().")
    end
    #---------------------------

end




"""
`StrongField.exportData( comps::Array{StrongField.Computation,1}, results::Array{Dict{String,Any},1}, dataLabel::String = "StrongFieldData" )`
    ... exports the results = [Array1 Array2 ...] returned by StrongField.perform with the StrongField computations comps = [computation1 computation2 ...] 
        into files with name dataLabel-dataType-1.csv, dataLabel-dataType-2.csv, ... where dataType = energy_distribution, azimuthal_angular_distribution, etc.
"""
function exportData( comps::Array{StrongField.Computation,1}, results::Array{Dict{String,Any},1}, dataLabel::String = "StrongFieldData" )
    for j = 1:size(results)[1]
        w = results[j]
        if typeof(comps[j].observable) == StrongField.SfaEnergyDistribution
            energyDistribution = w["energy distribution"]
            writedlm(dataLabel * "-energy_distribution" * "-" * string(j) * ".csv",hcat(energyDistribution.energies,energyDistribution.probabilities))
        elseif typeof(comps[j].observable) == StrongField.SfaMomentumDistribution
            angularDistribution = w["momentum distribution"]
            writedlm(dataLabel * "-momentum_distribution" * "-" * string(j) * ".csv",hcat(angularDistribution.phis,angularDistribution.probabilities))
        elseif typeof(comps[j].observable) == StrongField.SfaAzimuthalAngularDistribution
            angularDistribution = w["angular distribution"]
            writedlm(dataLabel * "-azimuthal_angular_distribution" * "-" * string(j) * ".csv",hcat(angularDistribution.phis,angularDistribution.probabilities))
        elseif typeof(comps[j].observable) == StrongField.SfaPolarAngularDistribution
            angularDistribution = w["angular distribution"]
            writedlm(dataLabel * "-polar_angular_distribution" * "-" * string(j) * ".csv",hcat(angularDistribution.phis,angularDistribution.probabilities))
        end
    end
end






#------------------------TEST: Plot the radial wave functions (JAC + hydrogenic)------------------------------------
"""
`StrongField.exportRadialWavefunctions( comps::Array{StrongField.Computation,1}, dataLabel::String = "StrongFieldData", savePlot::Bool = false )`
    ... exports the radial wave functions (initial state) = [Array1 Array2 ...] that are used in the StrongField computations comps = [computation1 computation2 ...] 
        into files with name dataLabel-radial_wavefunction-1.csv, dataLabel-radial_wavefunction-2.csv, etc.
        - If savePlot == true, the wave functions are also plotted and exported into a single figure dataLabel-radial_wave_function.pdf
"""
function exportRadialWavefunctions( comps::Array{StrongField.Computation,1}, dataLabel::String = "StrongFieldData", savePlot::Bool = false )

        minIonizationPotential = Float64

        for k = 1:size(comps)[1]

            #Extract the initial orbital of the active electron from the many-electron comp.initialLevel and set quantum numbers
            initialOrbitals = comps[k].initialLevel.basis.orbitals
            
            #Find highest lying orbital (smallest ionization potential)
            defaultSubshell = [sh for (sh,or) in initialOrbitals][1] #This is not nice; must be a better way to simply get a default element from a Dict
            o = initialOrbitals[defaultSubshell]
            minIonizationPotential = abs(o.energy)
            for (subshell,orbital) in initialOrbitals
                if   abs(orbital.energy) < minIonizationPotential
                    o = orbital
                    minIonizationPotential = abs(orbital.energy)
                end
            end
            
            ls = LevelSymmetry(o.subshell)
            n = o.subshell.n;      l = Int((ls.J.num+1)/2);    j = ls.J.num/2;
            
            if  (sign((-1)^l) == -1 && ls.parity == plus::Parity) || (sign((-1)^l) == 1 && ls.parity == minus::Parity)
                l = l - 1
            end
            l = floor(Int,l)
            
            if  comps[k].settings.hydrogenic
                if comps[k].settings.hydrogenic1s
                    P = StrongField.HydrogenPnl( o.energy, 1, 0, rGrid.r )
                else
                    P = StrongField.HydrogenPnl( o.energy, n, l, rGrid.r )
                end
            else
                P = o.P
            end
            
            writedlm(dataLabel * "-initial_radial_wavefunction" * string(k) * ".csv",hcat(comps[k].grid.r,P))
            
            if savePlot
                if k == 1
                    Plots.plot(comps[k].grid.r, P, 
                                                title = "Radial wave functions",
                                                xlabel = "r (a.u.)",
                                                ylabel = "P(r)",
                                                markershape = :circle,
                                                gridlinewidth = 2,
                                                tickfontsize = 10,
                                                labelfontsize = 10,
                                                labelfontfamily = "Latin Modern Roman",
                                                titlefontfamily = "Latin Modern Roman"
                                )
                else
                    Plots.plot!(comps[k].grid.r, P)
                end
            
            end
            
            
        end
        
        if savePlot                                    
            savefig(dataLabel * "-initial_radial_wave_function.pdf")
        end
end


#if false
#    minIonizationPotential = 0.
    #Extract the initial orbital of the active electron from the many-electron comp.initialLevel and set quantum numbers
#    initialOrbitals = initialLevel.basis.orbitals
    
    #Find highest lying orbital (smallest ionization potential)
#    defaultSubshell = [sh for (sh,or) in initialOrbitals][1] #This is not nice; must be a better way to simply get a default element from a Dict
#    o = initialOrbitals[defaultSubshell]
#    minIonizationPotential = abs(o.energy)
#    for (subshell,orbital) in initialOrbitals
#        if   abs(orbital.energy) < minIonizationPotential
#            global o = orbital
#            global minIonizationPotential = abs(orbital.energy)
#        end
#    end
    
#    ls = LevelSymmetry(o.subshell)
#    n = o.subshell.n;      l = Int((ls.J.num+1)/2);    j = ls.J.num/2;
#    
#    if  (sign((-1)^l) == -1 && ls.parity == plus::Parity) || (sign((-1)^l) == 1 && ls.parity == minus::Parity)
#        l = l - 1
#    end
#    l = floor(Int,l)
#    
#    hydrogenP = StrongField.HydrogenPnl( o.energy, n, l, rGrid.r )
#    
#    Plots.plot(rGrid.r, [o.P hydrogenP], 
#                                title = "Radial wave function: n=" * string(n) * ", l=" * string(l) * ", Ip=" * string(round(convertUnits("energy: from atomic to eV", o.energy),digits=2)) * " eV",
#                                xlabel = "r (a.u.)",
#                                ylabel = "P(r)",
#                                markershape = :circle,
#                                gridlinewidth = 2,
#                                tickfontsize = 10,
#                                labelfontsize = 10,
#                                labelfontfamily = "Latin Modern Roman",
#                                titlefontfamily = "Latin Modern Roman",
#                                label = ["JAC" "Hydrogen"]
#                              )
#    
#    savefig("radial_wavefunction.pdf")
#    
#    writedlm(dataName * "-radial_wavefunction.csv",hcat(rGrid.r,o.P,hydrogenP))
#end
##-----------------------------------------------------------------------------------------------------------------------


#if false
#epsilonp=2*omega
#kappa =-1
#lp=0

#nrContinuum = Continuum.gridConsistency(epsilonp, rGrid)
#contSettings = Continuum.Settings(false, nrContinuum)
#contSettings = Continuum.Settings(false, rGrid.NoPoints)

#newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, initialLevel.basis.subshells)
#newfLevel = Basics.generateLevelWithSymmetryReducedBasis(finalLevel, newiLevel.basis.subshells)
#newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, kappa), newiLevel)
#cOrbital, phase  = Continuum.generateOrbitalForLevel(epsilonp, Subshell(101, kappa), newfLevel, nuclearModel, rGrid, contSettings)

    
#Plots.plot( rGrid.r[1:size(cOrbital.P)[1]], real(cOrbital.P * exp(im*phase)), 
#                            title = "Radial wave function continuum",
#                            xlabel = "r (a.u.)",
#                            ylabel = "P(r)",
#                            markershape = :circle,
#                            gridlinewidth = 2,
#                            tickfontsize = 10,
#                            labelfontsize = 10,
#                            labelfontfamily = "Latin Modern Roman",
#                            titlefontfamily = "Latin Modern Roman"
#          )
#    
#    Z=1.0
#    cVolkov = StrongField.CoulombVolkovP( epsilonp, lp, Z, rGrid.r )
#    Plots.plot!( rGrid.r, real(cVolkov) )
#    
#    savefig("radial_wavefunction_continuum.pdf")


#end



