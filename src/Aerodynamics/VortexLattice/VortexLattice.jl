module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations

## Panel geometry and math tools
#==========================================================================================#

import ..AeroMDAO: Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform, p1, p2, p3, p4, dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient, weighted_vector, reshape_array, cartesian_to_freestream, freestream_to_cartesian

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")

export Horseshoe, horseshoe_line, horseshoe_point, bound_leg, bound_leg_center, bound_leg_vector, r1, r2, points, make_horseshoes

## Reference frames
#==========================================================================================#

include("reference_frames.jl")

export body_to_stability_axes, stability_to_body_axes, body_to_wind_axes, wind_to_body_axes, reflect_xz

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

export solve_system

"""
    make_horseshoes(horseshoe_panels)

Get the bound legs and collocation points of the horseshoes defined by horseshoe `Panel3D`s.
"""
make_horseshoes(horseshoe_panels) =	@. horseshoe_line(horseshoe_panels), horseshoe_point(horseshoe_panels)

"""
    solve_system(horseshoes, collocation_points, normals, U, Ω) 

Evaluate and return the vortex strengths ``\\Gamma``s given `Horseshoes`, their associated normal vectors (not necessarily the same as the panels' normals), the speed ``U`` and rotation vector ``\\Omega``.
"""
function solve_system(horseshoes, collocation_points, normals, U, Ω)
    V = map(r -> U + Ω × r, collocation_points)
    AIC  = influence_matrix(horseshoes, collocation_points, normals, -normalize(U))
    boco = boundary_condition(V, normals)
    Γs 	 = AIC \ boco 
end

## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics

## Post-processing
#==========================================================================================#

export case_dynamics, evaluate_coefficients

"""
    case_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

Compute the nearfield and farfield forces over a given component with associated vortex strengths `Γ_comp` and horseshoes `hs_comp`, using the vortex strengths ``\\Gamma``s and `horseshoes` of the entire configuration, the speed ``U``, angles of attack ``\\alpha`` and sideslip ``\\beta``, the rotation vector ``\\Omega``, the freestream density ``\\rho`` and the reference location for moments ``r_\\text{ref}``.
"""
function case_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)
    # Compute near-field dynamics
    surface_forces, surface_moments = nearfield_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, r_ref)

    # Compute farfield dynamics
    trefftz_force = trefftz_forces(Γ_comp, hs_comp, norm(U), α, β, ρ)

    surface_forces, surface_moments, trefftz_force
end

# In case of only one component
case_dynamics(Γs, horseshoes, U, α, β, Ω, ρ, r_ref) = case_dynamics(Γs, horseshoes, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

function evaluate_coefficients(forces, moments, trefftz_force, U, α, β, Ω, ρ, S, c, b)
    V = norm(U)
    q = dynamic_pressure(ρ, V)

    # Computing summed coefficients
    force, moment = sum(forces), sum(moments)

    # Transform near-field dynamics to wind axes
    trans_force	 = body_to_wind_axes(force, α, β)
    trans_force  = [ nearfield_drag(force, U); trans_force[2:end] ]
    trans_moment = body_to_wind_axes(stability_flip(moment), α, β)

    # Compute coefficients
    nearfield_coeffs = aerodynamic_coefficients(trans_force, trans_moment, Ω, V, S, b, c, ρ)
    farfield_coeffs	 = force_coefficient(trefftz_force, q, S)

    # Non-dimensional panel coefficients
    CFs = force_coefficient.(forces, q, S)
    CMs = moment_coefficient.(moments, q, S, b, c)

    nearfield_coeffs, farfield_coeffs, CFs, CMs
end

# Streamlines
include("streamlines.jl")

export streamlines

# EXPERIMENTAL: SYSTEM
include("system.jl")

export VLMState, VLMSystem, compute_influence_matrix!, compute_boundary_condition!, compute_collocation_points!, compute_horseshoes!, compute_normals!, solve_system!, compute_nearfield_forces!, compute_farfield_forces!, solve_case!, set_horseshoe_panels!, set_camber_panels!, compute_wake_properties!, normals, horseshoes, collocation_points, AIC, RHS, circulations, surface_force_coefficients, surface_moment_coefficients, total_force_coefficients, farfield_force_coefficients

end