# --- Bridson Poisson-disk sampler (3-D, periodic) ----------------------
# All coordinates here are plain `Float64`s in Å (no units) — units are
# stripped off at the boundary of `tip3p_water` and re-applied to the
# final `Atom` positions. Keeping internals unitless makes the rejection
# / distance-check inner loop straightforward and fast.

function _pbc_distance²(a::SVector{3, Float64},
                        b::SVector{3, Float64}, box::Float64)
   dx = a - b
   dx = SVector{3, Float64}(
      dx[1] - box * round(dx[1] / box),
      dx[2] - box * round(dx[2] / box),
      dx[3] - box * round(dx[3] / box),
   )
   return sum(abs2, dx)
end

function _far_from_all(c::SVector{3, Float64},
                       positions::Vector{SVector{3, Float64}},
                       d_min²::Float64, box::Float64)
   @inbounds for p in positions
      _pbc_distance²(c, p, box) < d_min² && return false
   end
   return true
end

# Phase 1: Bridson — k annulus tries per active point.
# Phase 2: uniform-rejection fallback to fill any holes Bridson missed.
# Throws if it can't reach `n_mol`.
function _poisson_disk_oxygens(box::Float64, n_mol::Int, d_min::Float64,
                               rng::AbstractRNG;
                               k::Int = 50,
                               uniform_fallback_trials::Int = 200_000)
   d_min²    = d_min * d_min
   positions = Vector{SVector{3, Float64}}(); sizehint!(positions, n_mol)
   active    = Int[]                         ; sizehint!(active,    n_mol)

   seed = SVector{3, Float64}(box * rand(rng), box * rand(rng), box * rand(rng))
   push!(positions, seed); push!(active, 1)

   # Phase 1: Bridson.
   while !isempty(active) && length(positions) < n_mol
      idx       = rand(rng, 1:length(active))
      center_id = active[idx]
      center    = positions[center_id]
      placed    = false
      for _ in 1:k
         u = SVector{3, Float64}(randn(rng), randn(rng), randn(rng))
         u = u / norm(u)
         r = d_min * (1.0 + rand(rng))             # uniform in [d_min, 2·d_min]
         c = center + r * u
         c = SVector{3, Float64}(mod(c[1], box), mod(c[2], box), mod(c[3], box))
         if _far_from_all(c, positions, d_min², box)
            push!(positions, c)
            push!(active, length(positions))
            placed = true
            break
         end
      end
      if !placed
         active[idx] = active[end]
         pop!(active)
      end
   end

   # Phase 2: uniform-rejection fallback.
   trials = 0
   while length(positions) < n_mol && trials < uniform_fallback_trials
      trials += 1
      c = SVector{3, Float64}(box * rand(rng), box * rand(rng), box * rand(rng))
      if _far_from_all(c, positions, d_min², box)
         push!(positions, c)
      end
   end

   length(positions) == n_mol ||
      error("Poisson-disk jammed: placed $(length(positions))/$n_mol oxygens " *
            "(box = $box Å, d_min = $d_min Å, fallback trials = $trials). " *
            "Lower d_min or n_mol.")
   return positions
end

# --- Uniform random rotation on SO(3) via the unit-quaternion method ---
function _random_rotation_matrix(rng::AbstractRNG)
   q = SVector{4, Float64}(randn(rng), randn(rng), randn(rng), randn(rng))
   q = q / norm(q)
   w, x, y, z = q[1], q[2], q[3], q[4]
   return SMatrix{3, 3, Float64}(
      1 - 2*(y*y + z*z),   2*(x*y + z*w),       2*(x*z - y*w),
      2*(x*y - z*w),       1 - 2*(x*x + z*z),   2*(y*z + x*w),
      2*(x*z + y*w),       2*(y*z - x*w),       1 - 2*(x*x + y*y),
   )
end
