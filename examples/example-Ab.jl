
println("Ab) Read in a full multiplet representation from files compatible to Ratip2012 and Grasp2018.")
#
wa = JAC.ManyElectron.Multiplet("from Ratip2012", "../test/files/ratip2012-at-is-a1-csl.inp", "../test/files/ratip2012-at-is-a1-scf.out", 
                                                  "../test/files/ratip2012-at-is-a1-relci.mix")

println("Not yet implemented; files do not exist !")
wa = JAC.ManyElectron.Multiplet("from Grasp2018", "../test/files/grasp2018-at-is-a1-csl.inp", "../test/files/grasp2018-at-is-a1-scf.out", 
                                                  "../test/files/grasp2018-at-is-a1-rci.mix")
