## Farfield dynamics
#==========================================================================================#

trefftz_potential(r_i, r_j, Γ_j) = let r = r_i - r_j; Γ_j / 2π * SVector(1, 0, 0) × r / norm(r)^2 end

∇φ(r, points, Γs) = sum(x -> trefftz_potential(r, x[1], x[2]), zip(points, Γs))

∂φ∂n(trefftz_line :: Line, points, Γs, normal) = dot(∇φ(center(trefftz_line), points, Γs), normal)

∂φ∂ns(trefftz_lines :: Vector{<: Line}, Δφs, normals) = ∂φ∂n.(trefftz_lines, (Ref ∘ points)(trefftz_lines), (Ref ∘ fwddiff)([ 0; -Δφs; 0 ]), normals)

function trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream)
	U_hat               = SVector(1, 0, 0)
	U_ref, free_ref     = Ref(U_hat), Ref(freestream)
	trefftz_lines       = @. body_to_wind_axes(bound_leg(horseshoes[end,:][:]), free_ref)

	trefftz_vectors     = @. vector(trefftz_lines)
	trefftz_proj_vecs   = @. trefftz_vectors - dot(U_ref, trefftz_vectors) * U_ref
	normals             = @. U_ref × trefftz_proj_vecs

	@timeit "Dihedrals" dihedrals = @. atan(getindex(trefftz_proj_vecs, 3), getindex(trefftz_proj_vecs, 2))
	@timeit "Projected Leg Norms" Δs = norm.(trefftz_proj_vecs)

	trefftz_lines, trefftz_proj_vecs, normals, dihedrals, Δs
end

"""
	trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, V, ρ, symmetry)

Compute the aerodynamic forces in the Trefftz plane given cumulative doublet strengths ``\\Delta \\phi``s, Trefftz panel lengths ``Δs``, doublet-normal directional derivatives ``\\partial \\phi / \\partial n``, Trefftz panel angles ``\\theta``s, the freestream speed and density ``V, \\rho``, and an option for symmetry.
"""
function trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, V, ρ, symmetry)     
	D_i = - 1/2 * ρ * sum(@. Δφs * Δs * ∂φ_∂n)
	Y   = - ρ * V * sum(@. Δφs * Δs * sin(dihedrals))
	L   = ρ * V * sum(@. Δφs * Δs * cos(dihedrals))
	
	symmetry ? SVector(D_i, 0, 2L) : SVector(D_i, Y, L)
end

"""
	trefftz_forces(Γs, horseshoes, freestream, ρ)

Compute the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream, ρ, symmetry :: Bool)
	# Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes
	trefftz_lines, trefftz_proj_vecs, normals, dihedrals, Δs = trefftz_preprocessing(horseshoes, freestream)

	# Compute directional derivatives of doublets in the normal direction
	@timeit "Sum Γs" Δφs  = vec(sum(Γs, dims = 1))
	@timeit "∂φ/∂n" ∂φ_∂n = ∂φ∂ns(trefftz_lines, Δφs, normals)

	# Compute forces
	trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, freestream.V, ρ, symmetry)
end 

"""
	farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

Compute farfield forces and moments in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, a reference point for calculation of moments, and a density ρ.
"""
function farfield_dynamics(Γs :: AbstractArray{<: Real}, horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225, symmetry = false)
	@timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, freestream, ρ, symmetry)
	@timeit "Trefftz Moment" trefftz_moment = body_to_wind_axes(r_ref, freestream) × trefftz_force

	trefftz_force, trefftz_moment
end