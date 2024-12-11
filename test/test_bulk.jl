
using AtomsBuilder
using Test, AtomsBase, Unitful, Random, JSON

##

@testset "Testing `bulk` and `repeat` against JuLIP reference data" begin
    test_systems = JSON.parsefile(joinpath(@__DIR__(), "..", "data", "test_systems.json"))

    _ustripvecvec(X) = [ ustrip.(x) for x in X ]

    compare_system(sys_f, sys_j) = (
          all( ustrip.( hcat(cell_vectors(sys_f)...) ) .≈ hcat(sys_j["cell"]...)' )  && 
          all( AtomsBuilder._convert_pbc(sys_j["pbc"]) .== periodicity(sys_f) ) && 
          all( atomic_number(sys_f, :) .== sys_j["Z"] ) && 
          all( _ustripvecvec(position(sys_f, :)) .≈ sys_j["X"] )    )

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
               sys_f = sys_f * nn
            end
            @test compare_system(sys_f, D["sys"])
        end
    end
end

@testset "Test argument errors are raised" begin
    @test_throws ArgumentError bulk(:Og)  # Case where no bulk structure is known
    @test_throws ArgumentError bulk(:Og)  # Case where no bulk structure is known
    @test_throws ArgumentError bulk(:Ti; cubic=true)  # Can't make hcp crystals cubic
end

@testset "Execute code" begin
    # not sure how to write a test for this, but at least it should execute
    sys0 = rattle!(bulk(:C, cubic=true) * (2,3,4), 0.1u"Å")
    sys1 = rattle!(bulk(:C, cubic=true) * (2,3,4), 0.1)
    sys2 = rattle!(bulk(:C) * (2,3,4), 0.1)
    rattle!(bulk(:C) * (2,3,4), 0.01u"nm")

    X = position(sys1, :)
    Xnew = [ x .+ 0.01u"Å" for x in X ]
    sys3 = set_positions(sys1, Xnew)
    @test position(sys3, :) == Xnew

    Z = atomic_number(sys1, :)
    @test all(Z .== AtomsBuilder.Chemistry.atomic_number(:C))
    zO = AtomsBuilder.Chemistry.atomic_number(:O)
    Znew = copy(Z); Znew[3:5:end] .= zO
    sys4 = set_elements(sys3, Znew)
    @test all(atomic_number(sys4, :) .== Znew)

    pold = deepcopy(sys4.particles) 
    deleteat!(sys4, 1:5)
    @test position.(pold[6:end]) == position(sys4, :)
    @test mass.(pold[6:end]) == mass(sys4, :)
    @test atomic_number.(pold[6:end]) == atomic_number(sys4, :)
end
