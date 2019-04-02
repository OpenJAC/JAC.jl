using Test
using JAC

@testset "Name" begin
    ## @test if 1 == 1  println("hallo");   global wa = true  else   println("Not hallo")   end

    global short = true
    printstyled("\nPerform tests on the JAC program; this may take a while ... \n", color=:cyan)
    JAC.define("print test: open", pwd() * "/runtests.report")

    include("test-halfintegers.jl")
    include("test-wignersymbols.jl")
    include("test-racahsum.jl")

    @test JAC.testMethods(short=short)
    @test JAC.testAmplitudes(short=short)
    @test JAC.testProperties(short=short)
    @test JAC.testProcesses(short=short)

    JAC.define("print test: close", "")
end
