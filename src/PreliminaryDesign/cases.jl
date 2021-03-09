"""
	solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, freestream :: Freestream, r_ref = [0.25, 0., 0.], ρ = 1.225; symmetry = false)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, with an option for symmetry.
"""
function solve_case(horseshoe_panels :: Matrix{<: Panel3D}, normals, freestream :: Freestream, ρ, r_ref = [0.25, 0., 0.]; symmetry = false)
	# Solve system
	Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], normals[:], freestream, symmetry), size(horseshoe_panels)...)

	# Compute forces and moments
	trans_forces, trans_moments, trans_rates, trefftz_force, trefftz_moment = case_dynamics(Γs, horseshoes, freestream, r_ref, ρ, symmetry)
	
	trans_forces, trans_moments, trans_rates, trefftz_force, trefftz_moment, Γs, horseshoes
end

"""
	solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref = [0.25, 0., 0.]; span_num :: Integer = 5, chord_num :: Integer = 10)

Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
"""
function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref = [0.25, 0., 0.]; span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer)
	# Determine spanwise panel distribution
	span_nums = ifelse(typeof(span_num) <: Integer, ceil.(Int, span_num / 2 .* wing.right.spans / span(wing.right)), span_num)

	# Compute panels and normals
	horseshoe_panels, camber_panels = vlmesh_wing(wing, span_nums, chord_num)
	normals = panel_normal.(camber_panels)

	# Compute forces and moments
	panel_forces, panel_moments, wind_rates, trefftz_force, trefftz_moment, Γs, horseshoes = solve_case(horseshoe_panels, normals, freestream, ρ, r_ref)

	# Compute aerodynamic coefficients
	nearfield_coeffs, farfield_coeffs, CFs, CMs  = case_coefficients(wing, panel_forces, panel_moments, wind_rates, trefftz_force, trefftz_moment, freestream.V, ρ)
	
	nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs
end

"""
	solve_case(foil :: Foil, uniform :: Uniform2D; sources = false, wake_length = 1e3, num_panels :: Integer = 60)

Evaluate a doublet-source case given a `Foil` with a `Uniform2D`, with optional named arguments to specify whether the source terms are non-zero, the length of the wake, and the number of panels for the analysis.
"""
function solve_case(foil :: Foil, freestream :: Uniform2D; viscous = true, sources = false, wake_length = 1e3, num_panels :: Integer = 60, wake_panels = 15)
	panels = paneller(foil, num_panels)
	cps, cls, cl_wake = solve_problem(panels, velocity(freestream), ifelse(viscous, wake_panels, sources), wake_length)
	# cdp, cl = sincos(freestream.ang) .* cl_wake
	# cm      = cl * r_ref
	# coeffs    = (cdp, cl, cm)

	cl_wake, cls, cps, panels
end

function case_coefficients(wing :: Union{Wing, HalfWing}, forces, moments, rates, trefftz_force, trefftz_moment, V, ρ)
	# Non-dimensionalisation quantities
	S = projected_area(wing)
	c = mean_aerodynamic_chord(wing)
	b = span(wing)

	# Non-dimensional dynamic coefficients
	CFs = force_coefficient.(forces, dynamic_pressure(ρ, V), S)
	CMs = moment_coefficient.(moments, dynamic_pressure(ρ, V), S, b, c)

	# Computing summed coefficients
	force, moment	 = sum(forces), sum(moments)
	nearfield_coeffs = aerodynamic_coefficients(force, moment, rates, V, S, b, c, ρ)
	farfield_coeffs	 = aerodynamic_coefficients(trefftz_force, trefftz_moment, rates, V, S, b, c, ρ)

	nearfield_coeffs, farfield_coeffs, CFs, CMs
end

symmetric_case_coefficients(wing :: Union{Wing, HalfWing}, force, moment, trans_rates, trefftz_force, trefftz_moment, V, ρ) = case_coefficients(wing, force, moment, trans_rates, trefftz_force, trefftz_moment, V, 2ρ)

# SYMMETRIC CASE WRONG RESULTS???
# if typeof(wing) == Wing && wing.left === wing.right && freestream.β == 0. && freestream.Ω == zeros(3)
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

# function solve_case(components :: Dict{<: Aircraft, Tuple{Vector{<: Integer}, Integer}}, freestream :: Freestream, ρ, r_ref = [0.25, 0., 0.])
#     # Compute panels
#     meshes = [ vlmesh_wing(component, panel_nums...) for (component, panel_nums) in components ] 
#     horseshoe_panels, camber_panels = first.(meshes), last.(meshes)
# 	normals = panel_normal.(camber_panels)

#     # Solve system
#    	Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], normals[:], freestream, symmetry), size(horseshoe_panels)...)

# 	# Compute forces on individual components
# end

# end
