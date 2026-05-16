
using AtomsBuilder
using AtomsBuilder.Examples
using Test, AtomsBase, Unitful, UnitfulAtomic, Random
using StaticArrays: SVector
using LinearAlgebra: norm

##

@testset "rocksalt — size, cell, pbc, charges" begin
   sys = nacl(2)

   @test length(sys) == 8 * 2^3
   @test periodicity(sys) == (true, true, true)

   cv = cell_vectors(sys)
   @test length(cv) == 3
   # Cubic box of side 2 · a_lat = 8 Å (default a_lat = 4 Å).
   for (k, v) in enumerate(cv)
      @test ustrip(u"Å", v[k]) ≈ 8.0
      for (j, vj) in enumerate(v)
         j == k && continue
         @test ustrip(u"Å", vj) ≈ 0.0 atol = 1e-12
      end
   end

   # All 64 atoms carry a `:charge` per-atom datum, ±1 e_au, summing to zero.
   charges = [a[:charge] for a in sys]
   @test all(c -> c == +1.0u"e_au" || c == -1.0u"e_au", charges)
   @test sum(charges) ≈ 0.0u"e_au"
   @test count(==( +1.0u"e_au"), charges) == 32
   @test count(==( -1.0u"e_au"), charges) == 32

   # Species: 4×4 of each per conventional cell, ×8 cells = 32 each.
   zs = atomic_number(sys, :)
   z_Na = AtomsBuilder.Chemistry.atomic_number(:Na)
   z_Cl = AtomsBuilder.Chemistry.atomic_number(:Cl)
   @test count(==(z_Na), zs) == 32
   @test count(==(z_Cl), zs) == 32
end

@testset "rocksalt — generic species + custom charge_label" begin
   sys = rocksalt(:Li, :F, 1; a_lat = 4.02u"Å", charge_label = :q)

   @test length(sys) == 8
   # Custom kwarg name flows through.
   for a in sys
      @test a[:q] == +1.0u"e_au" || a[:q] == -1.0u"e_au"
   end
   @test sum(a[:q] for a in sys) ≈ 0.0u"e_au"

   # Cubic box of side 1 · 4.02 Å.
   cv = cell_vectors(sys)
   @test ustrip(u"Å", cv[1][1]) ≈ 4.02
end

@testset "rocksalt — σ = 0 is deterministic" begin
   sys_a = nacl(2; rng = MersenneTwister(0))
   sys_b = nacl(2; rng = MersenneTwister(123))
   # No randomness when σ = 0 ⇒ positions identical regardless of seed.
   @test all(position(sys_a, :) .== position(sys_b, :))
end

@testset "rocksalt — σ > 0 displaces but preserves count and charges" begin
   sys0 = nacl(2)
   sys  = nacl(2; σ = 0.1u"Å", rng = MersenneTwister(42))

   @test length(sys) == length(sys0)
   # All atoms move by more than 0 (with probability 1 for a Gaussian).
   for i in 1:length(sys)
      @test norm(position(sys, i) - position(sys0, i)) > 0.0u"Å"
   end
   # Charges unchanged.
   @test sum(a[:charge] for a in sys) ≈ 0.0u"e_au"
end

@testset "rocksalt — rng reproducibility with σ > 0" begin
   sys_a = nacl(2; σ = 0.1, rng = MersenneTwister(7))
   sys_b = nacl(2; σ = 0.1, rng = MersenneTwister(7))
   @test all(position(sys_a, :) .== position(sys_b, :))

   sys_c = nacl(2; σ = 0.1, rng = MersenneTwister(8))
   @test position(sys_a, 1) != position(sys_c, 1)
end

@testset "rocksalt — bare-Real kwargs interpreted as Å" begin
   sys_q = nacl(1; a_lat = 4.0u"Å")
   sys_r = nacl(1; a_lat = 4.0)
   @test position(sys_q, :) == position(sys_r, :)
   @test cell_vectors(sys_q) == cell_vectors(sys_r)
end
