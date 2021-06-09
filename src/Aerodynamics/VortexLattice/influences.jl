# Matrix setup
#==========================================================================================#

"""
	influence_coefficient(r, horseshoe, normal, V_hat, symmetry)

Compute the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
influence_coefficient(horseshoe :: Horseshoe, r, normal, V_hat) = dot(velocity(r, horseshoe, 1., V_hat), normal)

"""
	influence_matrix(horseshoes, collocation_points, normals, V_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, and a unit vector representing the freestream.
"""
influence_matrix(horseshoes :: Vector{<: Horseshoe}, collocation_points, normals, V_hat) = [ influence_coefficient(horsie_j, r_i, n_i, V_hat) for (r_i, n_i) ∈ zip(collocation_points, normals), horsie_j ∈ horseshoes ]

"""
	boundary_condition(velocities, normals)

Compute the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)


## Pre-allocated versions
#====================================================#

function matrix_assembly!(AIC, RHS, horseshoes :: Vector{<: Horseshoe}, collocation_points, normals, V, Ω)
    for i ∈ 1:length(collocation_points)
        for j ∈ 1:length(horseshoes) 
            AIC[i,j] = dot(total_horseshoe_velocity(r1(collocation_points[i], horseshoes[j]), r2(collocation_points[i], horseshoes[j]), 1., -normalize(V)), normals[i])
        end
        RHS[i] = dot(V + Ω × collocation_points[i], normals[i])
    end
    nothing
end
    