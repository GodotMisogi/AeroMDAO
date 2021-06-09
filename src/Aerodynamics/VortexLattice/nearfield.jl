## Nearfield dynamics
#==========================================================================================#

kutta_joukowsky(ρ, Γ, V, l) = ρ * Γ * V × l
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r1(r, horseshoe), r2(r, horseshoe), Γ, V)

"""
    midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluate the induced velocity by the trailing legs at a given location ``r``, by summing over the velocities of Horseshoes with vortex strengths ``\\Gamma``s, rotation rates ``\\Omega``, and a freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r, horseshoes, Γs, U, Ω) = sum(x -> trailing_velocity(r, x[1], x[2], -normalize(U)), zip(horseshoes, Γs)) - U - Ω × r

nearfield_drag(force, U) = -dot(force, normalize(U))
horseshoe_moment(horseshoe :: Horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force
moments(horseshoes, forces, r_ref) = horseshoe_moment.(horseshoes, forces, Ref(r_ref))

"""
    nearfield_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ)

Compute the nearfield forces via the local Kutta-Jowkowski theorem given an array of horseshoes `hs_comp` to compute the forces on a component, their associated vortex strengths `Γ_comp`, the arrays of horseshoes and vortex strengths `Γs`  of the entire aircraft, the freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe, excluding the contribution of the bound leg.
"""
nearfield_forces(Γ_focus, hs_focus, Γs, horseshoes, U, Ω, ρ) = kutta_joukowsky.(Ref(ρ), Γ_focus, midpoint_velocity.(bound_leg_center.(hs_focus), Ref(horseshoes), Ref(Γs), Ref(U), Ref(Ω)), bound_leg_vector.(hs_focus))

function nearfield_dynamics(Γ_focus, hs_focus, Γs, horseshoes, U, Ω, ρ, r_ref)
    geom_forces  = nearfield_forces(Γ_focus, hs_focus, Γs, horseshoes, U, Ω, ρ)
    geom_moments = moments(hs_focus, geom_forces, r_ref)

    geom_forces, geom_moments
end