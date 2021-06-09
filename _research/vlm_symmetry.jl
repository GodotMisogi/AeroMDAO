
# SYMMETRIC CASE WRONG RESULTS???
# symmetric_case_coefficients(wing :: Union{Wing, HalfWing}, force, moment, trans_rates, trefftz_force, trefftz_moment, V, ρ) = case_coefficients(wing, force, moment, trans_rates, trefftz_force, trefftz_moment, V, 2ρ)

# if typeof(wing) == Wing && wing.left === wing.right && freestream.beta == 0. && freestream.omega == zeros(3)
# 	# Compute panels and normals
# 	horseshoe_panels, camber_panels = vlmesh_wing(wing.right, span_num, chord_num)
# 	normals = panel_normal.(camber_panels)

# 	panel_forces, panel_moments, wind_rates, force, moment, trefftz_force, trefftz_moment, Γs, horseshoes = solve_case(horseshoe_panels, normals, freestream, r_ref, ρ, symmetry = true)

# 	nearfield_coeffs, farfield_coeffs = symmetric_case_coefficients(wing.right, force, moment, wind_rates, trefftz_force, trefftz_moment, freestream.V, ρ)

# 	Γs					= reflect_mapper(identity, Γs)
# 	horseshoe_panels 	= reflect_mapper(x -> reflect_xz.(x), horseshoe_panels)
# 	camber_panels 		= reflect_mapper(x -> reflect_xz.(x), camber_panels)
# 	horseshoes 			= reflect_mapper(x -> reflect_xz.(x), horseshoes)
# else
# end

# if symmetry
#     col_vel = velocity(r, horseshoe, 1., V_hat)
#     ref_vel = (reflect_xz ∘ velocity)(reflect_xz(r), horseshoe, 1., V_hat)

#     dot(col_vel + ref_vel, normal)
# else

# if symmetry
# 	reflect_hs = reflect_xz.(horseshoes)
# 	geom_forces = [ nearfield_forces(Γs, reflect_hs, U, Ω, ρ)[end:-1:1]; 
# 					nearfield_forces(Γs, horseshoes, U, Ω, ρ) ]
# 	geom_moments = moments([ reflect_hs[end:-1:1]; horseshoes ], geom_forces, r_ref)
# else 