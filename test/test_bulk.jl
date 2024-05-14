
using AtomsBuilder, Test, AtomsBase
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

# not sure how to write a test for this: 
rattle!(bulk(:C, cubic=true) * (2,3,4), 0.1)