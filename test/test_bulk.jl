
using AtomsBuilder
using Test, AtomsBase, Unitful, Random, JSON
using AtomsSystems
using StaticArrays

##

@testset "Testing `bulk` and `repeat` against JuLIP reference data" begin
    test_systems = JSON.parsefile(joinpath(@__DIR__(), "..", "data", "test_systems.json"))

    for D in test_systems
        @testset "Testing $(D["sym"])" begin
            if D["sym"] in ("Ti", )
                cubic = false
            else
                cubic = D["cubic"]
            end

            sys_f = bulk(Symbol(D["sym"]); cubic, pbc=D["pbc"])
            nn = D["nn"] isa Integer ? D["nn"] : tuple(D["nn"]...)
            if nn != 1
               sys_f = repeat(sys_f, nn)
            end

            sys_j = D["sys"]
            @test cell_matrix(sys_f) ≈ hcat(sys_j["cell"]...)' *u"Å"
            @test all( sys_j["pbc"] .== periodicity(sys_f) )
            @test all( atomic_number(sys_f, :) .== sys_j["Z"] )
            @test all( position(sys_f, :) .≈ [ SVector{3}(x)*u"Å" for x in sys_j["X"] ] )
        end
    end
end

@testset "Test argument errors are raised" begin
    @test_throws ArgumentError bulk(:Og)  # Case where no bulk structure is known
    @test_throws ArgumentError bulk(:Og)  # Case where no bulk structure is known
    @test_throws ArgumentError bulk(:Ti; cubic=true)  # Can't make hcp crystals cubic
end

