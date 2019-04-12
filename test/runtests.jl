using Test
using JAC

@testset "Name" begin
    ##x global short = true
    printstyled("\nPerform tests on the JAC program; this may take a while ... \n", color=:cyan)
    JAC.define("print test: open", pwd() * "/runtests.report")

    include("inc-halfintegers.jl")
    include("inc-wignersymbols.jl")
    include("inc-racahsum.jl")

    ##x redirect_stdout(streamDummy) do
    @testset "JAC methods" begin
        @test JAC.testMethod_Wigner_3j() 
    end

    @testset "JAC amplitudes" begin
        @test JAC.testModule_MultipoleMoment() 
        @test JAC.testModule_ParityNonConservation() 
    end
    ##x close(streamDummy)

    @testset "JAC properties" begin
        ## @test JAC.testModule_Einstein()       ## printout too accurate
        @test JAC.testModule_Hfs()   
        @test JAC.testModule_LandeZeeman() 
        @test JAC.testModule_IsotopeShift()   
        @test JAC.testModule_AlphaVariation() 
        @test JAC.testModule_FormFactor() 
        @test JAC.testModule_DecayYield() 
        @test JAC.testModule_GreenFunction()
        @test JAC.testModule_MultipolePolarizibility()
        @test JAC.testModule_PlasmaShift() 
    end

    @testset "JAC processes" begin
        ## @test JAC.testModule_Radiative()     ## printout too accurate 
        @test JAC.testModule_PhotoExcitation()
        @test JAC.testModule_PhotoIonization()
        @test JAC.testModule_PhotoRecombination()
        ## @test JAC.testModule_AutoIonization()  
        @test JAC.testModule_Dielectronic()  
        ## @test JAC.testModule_PhotoExcitationFluores() 
        ## @test JAC.testModule_PhotoExcitationAutoion() 
        ## @test JAC.testModule_RayleighCompton() 
        ## @test JAC.testModule_MultiPhotonDeExcitation() 
        ## @test JAC.testModule_CoulombExcitation() 
    end

    ##x @test JAC.testMethods(short=short)
    ##x @test JAC.testAmplitudes(short=short)
    ##x @test JAC.testProperties(short=short)
    ##x @test JAC.testProcesses(short=short)

    JAC.define("print test: close", "")
end
