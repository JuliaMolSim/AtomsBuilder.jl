
using AtomsBuilder, Test, AtomsBase, Unitful, Random
import JuLIP

##

@info("Testing `bulk` and `repeat` against JuLIP reference implementation")

function _compare_particle(x1, x2)
   return (x1.position ≈ x2.position) && (x1.atomic_symbol == x2.atomic_symbol)
end

function _compare_at(at, at_f)
   return ( at.boundary_conditions == at_f.boundary_conditions &&
            at.bounding_box == at_f.bounding_box &&
            all(_compare_particle.(at.particles, at_f.particles)) )
end

for sym in [:Si, :Ge, :W, :Ti]
   at = bulk(sym)
   at_f = FlexibleSystem(JuLIP.bulk(sym))
   @test _compare_at(at, at_f)

   at = bulk(sym, cubic=true) 
   at_f = FlexibleSystem(JuLIP.bulk(sym, cubic=true))
   @test _compare_at(at, at_f)

   for nn in (rand(1:3), rand(2:4), 2)
      at = bulk(sym) * nn 
      at_f = FlexibleSystem(JuLIP.bulk(sym) * nn)
      @test _compare_at(at, at_f)

      at = bulk(sym, cubic=true) * nn
      at_f = FlexibleSystem(JuLIP.bulk(sym, cubic=true) * nn)
      @test _compare_at(at, at_f)
   end
end

##

# not sure how to write a test for this, but at least it should execute
sys0 = rattle!(bulk(:C, cubic=true) * (2,3,4), 0.1u"Å")
sys1 = rattle!(bulk(:C, cubic=true) * (2,3,4), 0.1)
sys2 = rattle!(bulk(:C) * (2,3,4), 0.1)

X = position(sys1)
Xnew = [ x .+ 0.01u"Å" for x in X ]
sys3 = set_positions(sys1, Xnew)
@test position(sys3) == Xnew

Z = atomic_number(sys1)
@test all(Z .== AtomsBuilder.Chemistry.atomic_number(:C))
zO = AtomsBuilder.Chemistry.atomic_number(:O)
Znew = copy(Z); Znew[3:5:end] .= zO
sys4 = set_elements(sys3, Znew)
@test all(atomic_number(sys4) .== Znew)

pold = deepcopy(sys4.particles) 
deleteat!(sys4, 1:5)
@test position.(pold[6:end]) == position(sys4)
@test atomic_mass.(pold[6:end]) == atomic_mass(sys4)
@test atomic_number.(pold[6:end]) == atomic_number(sys4)





