
println("Aa) Test of the radial integration in CI computations (Randolf).")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")

function printList(list)
  for element in list
    println(element)
  end
end

function getOrbital(orbitalList, subshell)
  for orbital in orbitalList
    if(subshell.n == orbital.subshell.n && subshell.kappa == orbital.subshell.kappa)
      return orbital
    end
  end
  
  throw(ErrorException("Orbital not found"))
  
end

nucleus = Nuclear.Model(10.)

grid = JAC.Radial.Grid(true)
setDefaults("standard grid", grid)


# configurations = [Configuration("1s^2 2s^2 2p^6")];      orbitalsIn = JAC.read("orbital list: Grasp92", "Work/Ne-jac/Ne-0+-1-scf.out")
# configurations = [Configuration("1s^2 2s^2 2p^5 3d")];   orbitalsIn = JAC.read("orbital list: Grasp92", "Work/Ne-jac/Ne-0+-2-scf.out")
configurations = [Configuration("1s^2 2s^2 2p^5")];      orbitalsIn = JAC.read("orbital list: Grasp92", "Work/Ne-jac/Ne-0+-3-scf.out")
# configurations = [Configuration("1s^2 2s^2 2p^5 3p")];   orbitalsIn = JAC.read("orbital list: Grasp92", "Work/Ne-jac/Ne-0+-4-scf.out")

# configurations = [Configuration("3s^2 3p^4 3d")]
println("Start from the non-relativistic configurations:")
printList(configurations)


# orbitals1 = JAC.read("orbital list: Grasp92", "Ne-0+-scf.exp.out")


relativisticConfigurationsList = ConfigurationR[]
for  configuration in configurations
    wa = JAC.generate("configuration list: relativistic", configuration)
    append!( relativisticConfigurationsList, wa)
end

println("Have the relativistic configurations:")
printList(relativisticConfigurationsList)

subshellList = JAC.generate("subshells: ordered list for relativistic configurations", relativisticConfigurationsList)
println("The subshells are: $subshellList")

setDefaults("relativistic subshell list", subshellList)

# Generate the relativistic CSF's for the given subshell list
csfList = CsfR[]
for  configuration in relativisticConfigurationsList
    newCsfs = generate("CSF list: from single ConfigurationR", configuration, subshellList)
    append!( csfList, newCsfs)
end
##x for  i = 1:length(csfList)    println("perform-ac: ", csfList[i])    end

println("Have the following $(length(csfList)) CSFs")
printList(csfList)

# Determine the number of electrons and the list of coreSubshells
NoElectrons      = sum( csfList[1].occupation )
coreSubshellList = Subshell[]
for  k in 1:length(subshellList)
    mocc = JAC.subshell_2j(subshellList[k]) + 1;    is_filled = true
    for  csf in csfList
        if  csf.occupation[k] != mocc    is_filled = false;    break   end
    end
    if   is_filled    push!( coreSubshellList, subshellList[k])    end
end
println("The number of electron is $NoElectrons")
println("The core subshells are: $coreSubshellList")




orbitals = Dict{Subshell, Orbital}()
for  subshell in subshellList
  try
    oin = getOrbital(orbitalsIn, subshell)
    orbitals = Base.merge( orbitals, Dict( subshell => oin) ) 
  catch err
    println("Not all orbitals found in input file, missing: $subshell")
    exit()
  end
end

# println(orbitals)


basis = (true, NoElectrons, subshellList, csfList, coreSubshellList, orbitals)

multiplet = perform("computation: CI", basis, nucleus, grid, MultipletSettings())


