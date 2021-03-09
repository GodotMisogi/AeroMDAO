module DoubletSource

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using TimerOutputs

using ..AeroMDAO: Panel2D, Point2D, collocation_point, point1, point2, trans_panel, affine_2D, panel_length, panel_angle, panel_tangent, panel_normal, panel_dist

## Non-dimensionalization
#===========================================================================#

include("../Tools/NonDimensional.jl")
import .NonDimensional: pressure_coefficient

## Math tools
#===========================================================================#

include("../Tools/MathTools.jl")
using .MathTools: rotation, inverse_rotation, midpair_map

## Solutions to Laplace's equation
#===========================================================================#

include("../Tools/Laplace.jl")
import .Laplace: source_potential, doublet_potential


## Matrix helpers
#===========================================================================#

function doublet_influence(panel_j :: Panel2D, panel_i :: Panel2D)
	xp, yp = trans_panel(panel_j, panel_i)
	doublet_potential(1., xp, yp, 0., panel_length(panel_j))
end

function source_influence(panel_j :: Panel2D, panel_i :: Panel2D)
	xp, yp = trans_panel(panel_j, panel_i)
	source_potential(1., xp, yp, 0., panel_length(panel_j))
end

boundary_condition(panel_j :: Panel2D, panel_i :: Panel2D, u) = -source_influence(panel_j, panel_i) * dot(u, panel_normal(panel_j))

function wake_panel(panels, bound)
	lastx, lasty = (point2 ∘ last)(panels)
	Panel2D(SVector(lastx, lasty), SVector(bound * lastx, lasty))
end

function wake_panels(panels, bound, num)
	lastx, lasty = (point2 ∘ last)(panels)
	bounds = range(lastx, bound, length = num)
	@. Panel2D(SVector(bounds[1:end-1], lasty), SVector(bounds[2:end], lasty))
end

## Dynamics helpers
#===========================================================================#

panel_velocity(dφ, dr, u, α) = dφ / dr + dot(u, α)

lift_coefficient(cp, dist_colpoints, panel_angle) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength, speed) = 2. * wake_strength / speed

## Matrix assembly
#===========================================================================#

include("matrix_func.jl")
include("matrix_prealloc.jl")

export solve_problem

function solve_problem(panels :: Vector{<: Panel2D}, u, sources :: Bool, wake_length)
	φs			= solve_strengths(panels, u, sources; bound = wake_length)
	cps, cls	= aerodynamic_coefficients(panels, φs, u, sources)
	cl_wake 	= lift_coefficient(last(φs), norm(u))

	cps, cls, cl_wake
end

function solve_problem(panels :: Vector{<: Panel2D}, u, num_wake :: Integer, wake_length)
	wakes 		= wake_panels(panels, wake_length, num_wake)
	φs			= solve_strengths(panels, u, wakes; bound = wake_length)
	cps, cls	= aerodynamic_coefficients(panels, φs, u, true)
	cl_wake 	= lift_coefficient(last(φs), norm(u))

	cps, cls, cl_wake
end

# Pre-allocated versions
# @timeit "Solve System (Pre-allocated)" 
# φs = solve_strengths_prealloc(panels, u)
# @timeit "Lift Coefficient (Pre-allocated)" 
# cl = lift_coefficient_prealloc(panels, φs, u)

end