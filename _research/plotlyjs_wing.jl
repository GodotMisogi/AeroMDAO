using PlotlyJS

panel_splits(coords) = [ [ c[1] for c in panel ] for panel in coords ], [ [ c[2] for c in panel ] for panel in coords ], [ [ c[3] for c in panel ] for panel in coords ]

trace_coords(xs, ys, zs, color = :black) = [    scatter3d(
                                                            x = x,
                                                            y = y,
                                                            z = z,
                                                            mode = :lines, 
                                                            line = attr(color = color),
                                                            showlegend = false,
                                                        ) 
                                                for (x, y, z) in zip(xs, ys, zs) ]

trace_panels(panels :: Vector{<: Panel3D}) = trace_coords(panel_splits(panel_coords.(panels))...)

trace_surface(wing :: Union{HalfWing, Wing}, span_num = 5, chord_num = 30; rotation = [1 0 0; 0 1 0; 0 0 1], translation = SVector(0, 0, 0)) = trace_panels([ transform(panel, rotation, translation) for panel in mesh_wing(wing, span_num, chord_num)[:] ])

function trace_panels(panels :: Vector{<: Panel3D}, Γs :: Vector{<: Real})
    coords = panel_coords.(panels[:])

    min_Γ, max_Γ = extrema(Γs)
    Γ_range = -map(-, min_Γ, max_Γ)
    norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ];

    xs, ys, zs = panel_splits(coords)

    trace = [   mesh3d(
                        x = x,
                        y = y,
                        z = z,
                        intensity = fill(norm_Γ, length(x)),
                        text = norm_Γ,
                        showscale = false,
                    ) 
                for (x, y, z, norm_Γ) in zip(xs, ys, zs, norm_Γs) ]
end

trace_streamlines(freestream :: Freestream, points, horseshoes :: Vector{<: Horseshoe}, Γs :: Vector{<: Real}, length :: Real, num_steps :: Integer = 100) = trace_coords((panel_splits ∘ streamlines)(freestream, points, horseshoes, Γs, length, num_steps)..., :lightblue)

function plot_case(horseshoe_panels, camber_panels, Γs, horseshoes, freestream, seed, length = 2, num_steps = 100)
trace_horsies = trace_panels(horseshoe_panels[:])
trace_horses  = trace_panels(horseshoe_panels[:], Γs[:])
trace_cambers = trace_panels(camber_panels[:])
trace_streams = trace_streamlines(freestream, seed, horseshoes[:], Γs[:], length, num_steps)

layout = Layout(scene = attr(aspectratio=attr(x=1,y=1,z=1)))
                plot = PlotlyJS.plot(
                                    [
                                        (trace for trace in trace_horsies)...,
                                        (trace for trace in trace_horses)...,
                                        (trace for trace in trace_streams)...,
                                        (trace for trace in trace_cambers)...,
                                    ],
                                    layout)

layout, plot
end

layout, plt = plot_case(horseshoe_panels, camber_panels, Γs, horseshoes, uniform, seed, 5., 100);
plt