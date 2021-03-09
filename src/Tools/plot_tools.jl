using CoordinateTransformations
using Rotations

plot_panels(panels :: Vector{<: Panel3D}) = (tupvector ∘ panel_coords).(panels)

function plot_wing(wing :: Union{HalfWing, Wing}; rotation = [1 0 0; 0 1 0; 0 0 1], translation = [0,0,0]) 
	leading, trailing = wing_bounds(wing)
	[ tuple((Translation(translation) ∘ LinearMap(rotation))(coords)...) for coords in [ leading; trailing[end:-1:1]; [ first(leading) ] ] ]
end

plot_streams(freestream, points, horseshoes, Γs, length, num_steps) = tupvector.(streamlines(freestream, points, horseshoes, Γs, length, num_steps))

plot_surface(wing :: Union{HalfWing, Wing}, span_num = 5, chord_num = 30; rotation = one(RotMatrix{3, Float64}), translation = SVector(0, 0, 0)) = plot_panels([ transform(panel, rotation, translation) for panel in mesh_wing(wing, span_num, chord_num)[:] ])


## Doublet-source
#==========================================================================================#

## Plotting domain
# x_domain, y_domain = (-1, 2), (-1, 1)
# grid_size = 50
# x_dom, y_dom = linspace(x_domain..., grid_size), linspace(y_domain..., grid_size)
# grid = x_dom × y_dom

# vels, pots = grid_data(dub_src_panels, grid)
# cp = pressure_coefficient.(uniform.mag, vels);

# lower_panels, upper_panels = split_panels(dub_src_panels);

# ## Airfoil plot
# plot( (first ∘ collocation_point).(upper_panels), (last ∘ collocation_point).(upper_panels), 
#         label = "Upper", markershape = :circle,
#         xlabel = "x", ylabel = "C_p")
# plot!((first ∘ collocation_point).(lower_panels), (last ∘ collocation_point).(lower_panels),
#         label = "Lower", markershape = :circle,
#         xlabel = "x", ylabel = "C_p")

# ## Pressure coefficient
# plot( (first ∘ collocation_point).(upper_panels), :cp .<< upper_panels, 
#         label = "Upper", markershape = :circle, 
#         xlabel = "x", ylabel = "C_p")
# plot!((first ∘ collocation_point).(upper_panels), :cp .<< lower_panels, 
#         label = "Lower", markershape = :circle, yaxis = :flip)

# ## Control volume
# p1 = contour(x_dom, y_dom, cp, fill = true)
# plot(p1)
# plot!(first.(:start .<< panels), last.(:start .<< panels), 
#       color = "black", label = "Airfoil", aspect_ratio = :equal)