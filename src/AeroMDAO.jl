module AeroMDAO

#----------------------IMPORTS--------------------------------#
using StaticArrays
using Rotations
using LinearAlgebra

using TimerOutputs

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
import .MathTools: tupvector, fwdsum, fwddiv, weighted_vector, vectarray, slope, splitat, adj3, columns,
cosine_dist, cosine_interp, # For Foil.jl
span, structtolist, inverse_rotation, rotation, affine_2D,
Point2D, Point3D, x, y, z  # For DoubletSource.jl


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
using .NonDimensional

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics, reynolds_number

## Panels
#===========================================================================#

include("Geometry/PanelGeometry.jl")
using .PanelGeometry

export Panel, Panel2D, Point2D, collocation_point


## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry.jl")
# using .AircraftGeometry

## Vortex lattice
#==========================================================================================#

include("VortexLattice/VortexLattice.jl")
using .VortexLattice

export Horseshoe, Freestream, velocity, streamlines, solve_horseshoes, transform, panel_coords

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

include("cases.jl")

export solve_case

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_surface, plot_streamlines, trace_surface, trace_panels, trace_coords, trace_streamlines, panel_splits

end