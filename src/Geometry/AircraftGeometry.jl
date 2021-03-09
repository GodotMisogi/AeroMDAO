module AircraftGeometry

using Base.Math
using Base.Iterators
using DelimitedFiles
using StaticArrays
using CoordinateTransformations
using Rotations
using LinearAlgebra

using ..AeroMDAO: cosine_dist, cosine_interp, splitat, adj3, slope, columns, fwdsum, fwddiv, weighted_vector, vectarray, Point2D, Panel2D, Panel3D, extend_yz, transform

abstract type Aircraft end

export Aircraft

## Foil geometry
#==========================================================================================#

include("foil.jl")

export Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, cosine_foil, camthick_to_CST, coords_to_CST

## Wing geometry
#==========================================================================================#

include("wing.jl")

export HalfWing, Wing, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, print_info, leading_edge, leading_chopper, trailing_chopper, wing_chopper, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing

end