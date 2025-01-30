using AtomsBase
using AtomsBuilder
using Test

@testset "PubChem inteface" begin
    sys = load_from_pubchem( "water" )
    @test length(sys) == 3
    sys = load_from_pubchem( 887 ) # methanol
    @test length(sys) == 6
    sys = load_from_pubchem( smiles="CC(=O)C" ) # asetone
    @test length(sys) == 10
    sys = load_from_pubchem( "64-17-5" ) # ethanol
    @test length(sys) == 9
end