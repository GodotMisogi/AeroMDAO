module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs

## Math tools
#==========================================================================================#

include("../Tools/MathTools.jl")

using .MathTools: accumap, fwddiff, structtolist, quarter_point, three_quarter_point

## Freestream
#==========================================================================================#

include("../Tools/Laplace.jl")
import .Laplace: Freestream, velocity, aircraft_velocity

export Freestream, velocity, aircraft_velocity

## Panel geometry
#==========================================================================================#

import ..AeroMDAO: Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform, point1, point2, point3, point4

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")

export Horseshoe, horseshoe_line, horseshoe_point

## Reference frames
#==========================================================================================#

include("../Tools/ReferenceFrames.jl")
# using .ReferenceFrames: body_to_stability_axes, body_to_wind_axes

export body_to_stability_axes, body_to_wind_axes, reflect_xz

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

export solve_horseshoes

make_horseshoes(horseshoe_panels) =	@. horseshoe_line(horseshoe_panels), horseshoe_point(horseshoe_panels)

function solve_system(colpoints, horseshoes, normals, total_vel, U, symmetry)
	AIC = influence_matrix(colpoints, normals, horseshoes, -normalize(U), symmetry)
	boco = boundary_condition(total_vel, normals)

	AIC \ boco
end

"""
	solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry) 

Evaluate and return the vortex strengths ``\\Gamma s`` and associated horsehoes given Panel3Ds, their associated normal vectors (not necessarily the same as the panels' normals), and a Freestream, with the option to use the symmetry of the problem in the ``x``-``z`` plane.
"""
function solve_horseshoes(horseshoe_panels, normals, freestream :: Freestream, symmetry = false) 
	U = aircraft_velocity(freestream)
	horseshoes, colpoints = make_horseshoes(horseshoe_panels)
	total_vel = [ U + freestream.Ω × colpoint for colpoint in colpoints ]
	Γs = solve_system(colpoints, horseshoes, normals, total_vel, U, symmetry)

	Γs, horseshoes
end


## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics

## Post-processing
#==========================================================================================#

export case_dynamics

function case_dynamics(Γs, horseshoes, freestream :: Freestream, r_ref, ρ, symmetry = false)
	# Compute near-field dynamics
	trans_forces, trans_moments, trans_rates = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ, symmetry)

	# Compute farfield dynamics
	trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ, symmetry)

	trans_forces, trans_moments, trans_rates, trefftz_force, trefftz_moment
end

# Streamlines
include("streamlines.jl")

export streamlines

end