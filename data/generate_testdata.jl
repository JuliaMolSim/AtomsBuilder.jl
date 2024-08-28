
# this requires to be run in an environment that has JuLIP and JSON
# installed. 
using JuLIP, JSON 

make_dict(sym, cubic, pbc, nn, at) = 
   Dict("sym" => sym, "cubic" => cubic, "pbc" => pbc, "nn" => nn, 
        "sys" => write_dict(at))

test_systems = []    

for sym in [:Si, :Ge, :W, :Ti]
   at = bulk(sym)
   push!(test_systems, make_dict(sym, false, true, 1, at))

   at = bulk(sym, cubic=true) 
   push!(test_systems, make_dict(sym, true, true, 1, at))
end

for ntest = 1:30
   sym = rand([:Si, :Ge, :W, :Ti]) 
   cubic = rand(Bool)
   pbc = tuple(rand(Bool, 3)...)
   nn = (rand(1:3), rand(2:4), 2)

   at = bulk(sym, cubic=cubic, pbc = pbc) * nn 
   push!(test_systems, make_dict(sym, cubic, pbc, nn, at))
end

open(@__DIR__() * "/test_systems.json", "w") do f
   JSON.print(f, test_systems)
end

# read the file back 
# test_systems2 = JSON.parsefile("test_systems.json")