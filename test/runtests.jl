using SeqFold
using Test
using Aqua

@testset "SeqFold.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SeqFold; ambiguities = false,)
    end
    # Write your tests here.
end
