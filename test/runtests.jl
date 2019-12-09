using Test
using JAC, ..Defaults, ..TestFrames

@testset "Name" begin
    ##x global short = true
    printstyled("\nPerform tests on the JAC program; this may take a while ... \n", color=:cyan)
    ## Defaults.Constants.define("print test: open", pwd() * "/runtests.report")

    ##x include("inc-halfintegers.jl")
    ##x include("inc-wignersymbols.jl")
    ##x include("inc-racahsum.jl")

    ##x redirect_stdout(streamDummy) do
    @testset "JAC methods" begin
        @test TestFrames.testMethod_Wigner_3j() 
    end

    @testset "JAC amplitudes" begin
        @test TestFrames.testModule_MultipoleMoment() 
        @test TestFrames.testModule_ParityNonConservation() 
    end
    ##x close(streamDummy)

    @testset "JAC properties" begin
        ## @test TestFrames.testModule_Einstein()       ## printout too accurate
        @test TestFrames.testModule_Hfs()   
        @test TestFrames.testModule_LandeZeeman() 
        @test TestFrames.testModule_IsotopeShift()   
        @test TestFrames.testModule_AlphaVariation() 
        @test TestFrames.testModule_FormFactor() 
        # @test TestFrames.testModule_DecayYield()  ... need to be adapted
        ## @test TestFrames.testModule_GreenFunction()
        @test TestFrames.testModule_MultipolePolarizibility()
        @test TestFrames.testModule_PlasmaShift() 
    end

    @testset "JAC processes" begin
        ##  @test TestFrames.testModule_PhotoEmission()     ## printout too accurate 
        @test TestFrames.testModule_PhotoExcitation()
        ## @test TestFrames.testModule_PhotoIonization()
        ## @test TestFrames.testModule_PhotoRecombination()
        ## @test TestFrames.testModule_AutoIonization()  
        @test TestFrames.testModule_Dielectronic()  
        ## @test TestFrames.testModule_PhotoExcitationFluores() 
        ## @test TestFrames.testModule_PhotoExcitationAutoion() 
        ## @test TestFrames.testModule_RayleighCompton() 
        ## @test TestFrames.testModule_MultiPhotonDeExcitation() 
        ## @test TestFrames.testModule_CoulombExcitation() 
    end

    ## Defaults.setDefaults("print test: close", "")
end
