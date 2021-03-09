module AeroMDAO

#----------------------IMPORTS--------------------------------#
using StaticArrays
using Rotations
using LinearAlgebra

using TimerOutputs

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
import .MathTools: tupvector, fwdsum, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns, extend_yz, reflect_mapper, cosine_dist, cosine_interp, structtolist, inverse_rotation, rotation, affine_2D, Point2D, Point3D, x, y, z

export tupvector


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
using .NonDimensional

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics, reynolds_number

## Panels
#===========================================================================#

include("Geometry/PanelGeometry.jl")
using .PanelGeometry

export Panel, Panel2D, Panel3D, Point2D, collocation_point, point1, point2, point3, point4, transform, panel_normal, midpoint, panel_location, panel_tangent


## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry.jl")
using .AircraftGeometry

export Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, cosine_foil, camthick_to_CST, coords_to_CST, # 2D setups
HalfWing, Wing, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, print_info, leading_edge, leading_chopper, trailing_chopper, wing_chopper, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing

## Vortex lattice
#==========================================================================================#

include("VortexLattice/VortexLattice.jl")
using .VortexLattice

export Horseshoe, Freestream, velocity, streamlines, solve_horseshoes, transform, horseshoe_line, horseshoe_point

## Doublet-source panel method
#==========================================================================================#

include("Tools/Laplace.jl")
import .Laplace: Uniform2D, velocity

export Uniform2D, velocity

include("DoubletSource/DoubletSource.jl")
using .DoubletSource

export lift_coefficient

## Aerodynamic analyses
#==========================================================================================#

include("PreliminaryDesign/cases.jl")

export solve_case

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_streams, plot_wing, plot_surface, plot_streamlines, trace_surface, trace_panels, trace_coords, trace_streamlines, panel_splits, plot_case

end