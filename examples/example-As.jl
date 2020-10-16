
#
println("As) Test of spin-angular coefficients.")

configs = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]
        
# Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
relconfList = ConfigurationR[]
for  conf in configs    wa = Basics.generate("configuration list: relativistic", conf);     append!( relconfList, wa)   end
for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end
subshellList = Basics.generate("subshells: ordered list for relativistic configurations", relconfList)
Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)

# Generate the relativistic CSF's for the given subshell list
csfList = CsfR[]
for  relconf in relconfList
    newCsfs = Basics.generate("CSF list: from single ConfigurationR", relconf, subshellList);   append!( csfList, newCsfs)
end


if  true
    # Calculate angular coefficients for a scalar one- or two-particle operator
    op = SpinAngular.OneParticleOperator(0, plus, true)     ## SpinAngular.TwoParticleOperator(0, plus, true)
    for  leftCsf in csfList
        for rightCsf in csfList
            coeffs = SpinAngular.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
        end
    end
    #
elseif  false
    # Calculate angular coefficients for a nonscalar one- particle operator
    for  leftCsf in csfList
        for rightCsf in csfList
            op     = SpinAngular.OneParticleOperator(1, plus, true)  ## SpinAngular.OneParticleOperator(1, minus, true)
            coeffs = SpinAngular.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
        end
    end
    #
end
