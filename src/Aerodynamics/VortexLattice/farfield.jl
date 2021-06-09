## Farfield dynamics
#==========================================================================================#

"""
    trefftz_velocity(r_i, r_j, Γ_j)

Compute the induced velocity at point ``r_i`` of the wake in the Trefftz plane due to the strength ``\\Gamma_j`` at ``r_j``. 
"""
trefftz_velocity(r_i, r_j, Γ_j) = let r = r_i - r_j; Γ_j / 2π * SVector(1, 0, 0) × r / norm(r)^2 end

"""
    trefftz_influence_matrix(centers, normals, points)

Compute the aerodynamic influence coefficient matrix of the wake in the Trefftz plane given the center points, normal vectors, and the points of the wake panels.
"""
trefftz_influence_matrix(centers, normals, points) = [ dot(trefftz_velocity(r_i, r_j2, 1.) - trefftz_velocity(r_i, r_j1, 1.), n_i) for (r_i, n_i) in zip(centers, normals), (r_j2, r_j1) in zip(points[2:end], points) ]


"""
    doublet_normal_derivatives(wake_lines :: Vector{<: Line}, Δφs, normals)

Compute the normal derivative strengths of the doublets given the wake `Line`s, net doublet strengths ``\\Delta \\phis``, and associated normal vectors.
"""
function doublet_normal_derivatives(wake_lines :: Vector{<: Line}, Δφs, normals)
    centers = center.(wake_lines)
    pts 	= points(wake_lines)
    AIC   	= trefftz_influence_matrix(centers, normals, pts)
    ∂φ_∂n 	= AIC * Δφs
end

project_vector(vector, U) = vector - dot(U, vector) * U
normal(wake_proj_vector, U) = normalize(U × wake_proj_vector)
dihedral(wake_proj_vec)  = atan(wake_proj_vec[3], wake_proj_vec[2])

"""
    trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, α, β)

Project `Horseshoe`s into the Trefftz plane aligned with the wind axes angles ``\\alpha,~\\beta``, and compute normal vectors, projection angles and lengths.
"""
function trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, α, β)
    # Reference velocity for broadcasting
    U_ref = (Ref ∘ SVector)(1, 0, 0)

    # Transform to wind axes
    wake_lines     = @. body_to_wind_axes(bound_leg(horseshoes[end,:][:]), α, β)

    # Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes
    wake_vectors   = @. vector(wake_lines); 
    wake_proj_vecs = @. project_vector(wake_vectors, U_ref)

    # Normals, dihedral angles, and projection lengths
    wake_normals   = @. normal(wake_proj_vecs, U_ref)
    wake_dihedrals = @. dihedral(wake_proj_vecs)
    wake_lengths   = @. norm(wake_proj_vecs)

    wake_lines, wake_normals, wake_dihedrals, wake_lengths
end

"""
    trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, V, ρ)

Compute the aerodynamic forces in the Trefftz plane given cumulative doublet strengths ``\\Delta \\phi``s, Trefftz panel lengths ``Δs``, doublet-normal directional derivatives ``\\partial \\phi / \\partial n``, Trefftz panel angles ``\\theta``s, the freestream speed and density ``V, \\rho``.
"""
function trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, V, ρ)     
    D_i = - 1/2 * ρ * sum(@. Δφs * Δs * ∂φ_∂n)
    Y   = - ρ * V * sum(@. Δφs * Δs * sin(θs))
    L   = ρ * V * sum(@. Δφs * Δs * cos(θs))
    
    SVector(D_i, Y, L)
end

# ifelse(symmetry, SVector(D_i, 0, 2L), SVector(D_i, Y, L))

"""
    trefftz_forces(Γs, horseshoes, freestream, ρ)

Compute the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{<: Horseshoe}, speed, α, β, ρ)
    # Get projections of horseshoes into Trefftz plane with the associated normals, dihedral angles and lengths
    wake_lines, normals, dihedrals, Δs = trefftz_preprocessing(horseshoes, α, β)

    # Compute directional derivatives of doublets in the normal direction
    Δφs   = vec(sum(Γs, dims = 1))
    ∂φ_∂n = doublet_normal_derivatives(wake_lines, Δφs, normals)

    # Compute forces
    trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, speed, ρ)
end