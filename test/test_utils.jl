
using AtomsBuilder, Test, AtomsBase, Unitful, Random

##

Random.seed!(1234)
sys = repeat(bulk(:Ti, cubic=false), 5)
sys = randz!(sys, [ :Ti => 0.5, :O => 0.5 ])
Z = atomic_number(sys, :)
@test count( Z .== 8 ) / length(Z) > 0.4
@test count( Z .== 22 ) / length(Z) > 0.4

