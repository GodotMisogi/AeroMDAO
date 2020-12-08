## 
using Revise
using StaticArrays
using BenchmarkTools
using TimerOutputs
using ProfileView
using AeroMDAO

## Wing section setup
foil = naca4((4,4,1,2))
num_secs = 3
foils = [ foil for i ∈ 1:num_secs ]

airfoils = Foil.(foils)
wing_chords = [0.18, 0.16, 0.08]
wing_twists = [2., 0., -2.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0., 11.3]
wing_sweeps = [1.14, 8.]

# wing_right = HalfWing(airfoils, wing_chords, wing_twists, wing_spans, wing_dihedrals, wing_sweeps)
wing_right = HalfWing(Foil.(naca4((2,4,1,2)) for i ∈ 1:5),
                        [0.0639628599561049, 0.06200381820887121, 0.05653644812231768, 0.04311297779068357, 0.031463501535620116],
                      [0., 0., 0., 0., 0.],
                      [0.2, 0.2, 0.2, 0.2],
                      [0., 0., 0., 0.],
                      [0., 0., 0., 0.])
wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 0.0, 0.0, Ω)
@time coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 5, print = true) 

print_timer();

## Panel method: TO DO
# wing_panels = mesh_wing(wing, 10, 30);
# wing_coords = plot_panels(wing_panels[:])
camber_coords = plot_panels(camber_panels[:])
horseshoe_coords = plot_panels(horseshoe_panels[:]);

## Streamlines
reset_timer!()

@timeit "Computing Streamlines" streams = plot_streamlines.(streamlines(uniform, horseshoe_panels[:], horseshoes, Γs, 2, 100));

print_timer()

##
min_Γ, max_Γ = extrema(Γs)
Γ_range = -map(-, min_Γ, max_Γ)
norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ]

##
using PlotlyJS

##
horse_xs = [ [ c[1] for c in panel ] for panel in horseshoe_coords ]
horse_ys = [ [ c[2] for c in panel ] for panel in horseshoe_coords ]
horse_zs = [ [ c[3] for c in panel ] for panel in horseshoe_coords ]

camber_xs = [ [ c[1] for c in panel ] for panel in camber_coords ]
camber_ys = [ [ c[2] for c in panel ] for panel in camber_coords ]
camber_zs = [ [ c[3] for c in panel ] for panel in camber_coords ]

streams_xs = [ [ c[1] for c in panel ] for panel in streams ]
streams_ys = [ [ c[2] for c in panel ] for panel in streams ]
streams_zs = [ [ c[3] for c in panel ] for panel in streams ];

# wing_xs = [ [ c[1] for c in panel ] for panel in wing_coords ]
# wing_ys = [ [ c[2] for c in panel ] for panel in wing_coords ]
# wing_zs = [ [ c[3] for c in panel ] for panel in wing_coords ];

##
layout = Layout(
                title = "Penguins",
                scene=attr(aspectmode="manual", aspectratio=attr(x=1,y=1,z=1)),
                zlim=(-0.1, 5.0)
                )

trace_horses = [ mesh3d(
                        x = x,
                        y = y,
                        z = z,
                        intensity = repeat([norm_Γ], length(x)),
                        text = norm_Γ,
                        showscale = false,
                        ) for (x, y, z, norm_Γ) in zip(horse_xs, horse_ys, horse_zs, norm_Γs) ]

trace_horsies = [ PlotlyJS.scatter3d(
                            x = x,
                            y = y,
                            z = z,
                            mode = :lines, 
                            line = attr(color = :black),
                            showlegend = false,
                            ) for (x, y, z) in zip(horse_xs, horse_ys, horse_zs) ]

trace_cambers = [ PlotlyJS.scatter3d(
                       x = x,
                       y = y,
                       z = z,
                       mode = :lines, 
                       line = attr(color = :black),
                       showlegend = false,
                       ) for (x, y, z) in zip(camber_xs, camber_ys, camber_zs) ]

trace_streams = [ PlotlyJS.scatter3d(
                            x = x, 
                            y = y, 
                            z = z, 
                            mode = :lines, 
                            line = attr(color = :lightblue),
                            showlegend = false,
                            ) for (x, y, z) in zip(streams_xs, streams_ys, streams_zs) ];

# trace_wing =    [ PlotlyJS.scatter3d(
#                             x = x,
#                             y = y,
#                             z = z,
#                             mode = :lines, 
#                             line = attr(color = :black),
#                             showlegend = false,
#                             ) for (x, y, z) in zip(wing_xs, wing_ys, wing_zs) ];

PlotlyJS.plot([ 
        [ trace for trace in trace_horses ]...,
        [ trace for trace in trace_horsies ]..., 
        [ trace for trace in trace_cambers ]...,
        [ trace for trace in trace_streams ]...,
        # [ trace for trace in trace_wing ]...,
     ], 
     layout)

##
# using Plots
# plotlyjs()

# ##
# plot(xaxis = "x", yaxis = "y", zaxis = "z", aspectratio = 1., size=(1280, 720))
# plot!.(camber_coords, color = :black, label = :none)
# [ mesh3d!(coord, colorscale = :viridis) for (coord, norm_Γ) in zip(horseshoe_coords, norm_Γs) ]
# plot!.(streams, color = :green, label = :none)

# gui();