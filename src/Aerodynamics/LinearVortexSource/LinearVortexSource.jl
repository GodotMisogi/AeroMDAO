module LinearVortexSource

using LinearAlgebra
using StaticArrays
using Base.Iterators

using ..AeroMDAO: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, Point2D, collocation_point, p1, p2, transform_panel, affine_2D, panel_length, panel_angle, panel_tangent, panel_normal, panel_dist, rotation, inverse_rotation, midpair_map, panel_velocity, pressure_coefficient, wake_panel, wake_panels, panel_points

include("singularities.jl")

export constant_source_velocity, linear_source_velocity_a, linear_source_velocity_b, linear_vortex_velocity_a, linear_vortex_velocity_b

include("matrix.jl")

export total_velocity, source_velocity, vortex_velocity, vortex_influence_matrix, source_influence_matrix, neumann_boundary_condition, kutta, two_point_matrix, source_matrix, vortex_matrix, constant_source_matrix, constant_source_boundary_condition

end