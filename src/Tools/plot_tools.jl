function plot_panels(panels :: Vector{<: Panel3D})
    coords = panel_coords.(panels)
    tupvector.([coord; [coord[1]]] for coord in coords)
end

function plot_wing(wing :: Union{HalfWing, Wing}, rotation, translation)
    affine = Translation(translation) ∘ LinearMap(rotation)
    leading, trailing = wing_bounds(wing)
    wing_coords = [ leading; trailing[end:-1:1]; [ first(leading) ] ]
    # foil_coords = [ [ [coord[1]; 0; coord[2]] .* chord .+ loc for coord in foil.coords ] for (chord, foil, loc) in zip(wing.right.chords[end:-1:1], wing.right.foils[end:-1:1], wing_coords) ]

    [ tuple(affine(coords)...) for coords in wing_coords ]
end

plot_wing(wing :: Union{HalfWing, Wing}; angle = 0., axis = [1., 0., 0.], position = zeros(3)) = plot_wing(wing, AngleAxis{Float64}(angle, axis...), position)

plot_streams(freestream, points, horseshoes, Γs, length, num_steps) = tupvector.(streamlines(freestream, points, horseshoes, Γs, length, num_steps))

plot_surface(wing :: Union{HalfWing, Wing}, span_num = 5, chord_num = 30; rotation = one(RotMatrix{3, Float64}), translation = SVector(0, 0, 0)) = plot_panels(transform(panel, rotation, translation) for panel in mesh_wing(wing, span_num, chord_num)[:])


## Doublet-source
#==========================================================================================#

## Plotting domain
# x_domain, y_domain = (-1, 2), (-1, 1)
# grid_size = 50
# x_dom, y_dom = linspace(x_domain..., grid_size), linspace(y_domain..., grid_size)
# grid = x_dom × y_dom

# vels, pots = grid_data(dub_src_panels, grid)
# cp = pressure_coefficient.(uniform.magnitude, vels);

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