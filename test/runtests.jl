using AtomsBuilder
using Test

@testset "AtomsBuilder.jl" begin
    # Write your tests here.
    @testset "bulk" begin include("test_bulk.jl"); end
    @testset "utils" begin include("test_utils.jl"); end
    @testset "Examples — rocksalt" begin include("test_examples_rocksalt.jl"); end
    @testset "Examples — water"    begin include("test_examples_water.jl");    end
    include("test_pubchem.jl")
end


