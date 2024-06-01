using AtomsBuilder
using Test

@testset "AtomsBuilder.jl" begin
    # Write your tests here.
    @testset "bulk" begin include("test_bulk.jl"); end
    @testset "utils" begin include("test_utils.jl"); end
end


