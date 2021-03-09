# Matrix setup
#==========================================================================================#

"""
	influence_coefficient(r, horseshoe, normal, V_hat, symmetry)

Compute the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
function influence_coefficient(r, horseshoe :: Horseshoe, normal, V_hat, symmetry :: Bool)
	if symmetry
		col_vel = velocity(r, horseshoe, 1., V_hat)
		ref_vel = (reflect_xz ∘ velocity)(reflect_xz(r), horseshoe, 1., V_hat)

		dot(col_vel + ref_vel, normal)
	else
		dot(velocity(r, horseshoe, 1., V_hat), normal)
	end
end

"""
	influence_matrix(colpoints, normals, horseshoes, V_hat, symmetry)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, and a unit vector representing the freestream.
"""
influence_matrix(colpoints, normals, horseshoes :: Vector{<: Horseshoe}, V_hat, symmetry :: Bool) = [ influence_coefficient(r_i, horsie_j, n_i, V_hat, symmetry) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]

"""
	boundary_condition(velocities, normals)

Compute the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)

# Pre-allocation attempts
#==========================================================================================#

function matrix_assembly!(AIC :: AbstractMatrix{<: Real}, boco, colpoints :: AbstractArray{<: SVector{3, <: Real}}, normals :: AbstractArray{<: SVector{3, <: Real}}, horseshoes :: Vector{<: Horseshoe}, U, Ω, V_hat)
	for i ∈ 1:length(colpoints)
		for j ∈ 1:length(horseshoes)
			
			AIC[i,j] = dot(total_horseshoe_velocity(r1(colpoints[i], horseshoes[j]), r2(colpoints[i], horseshoes[j]), 1., V_hat), normals[i])
		end
		boco[i] = dot(velocities[i], normals[i])
	end
	nothing
end
