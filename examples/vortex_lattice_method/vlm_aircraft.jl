## Aircraft analysis case
using AeroMDAO

## Wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [2.0, 0.0, -2.0],
            spans     = [4.0, 0.2],
            dihedrals = [5., 30.],
            sweep_LEs = [5., 30.]);
print_info(wing, "Wing")

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39])
print_info(htail, "Horizontal Tail")

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97])
print_info(vtail, "Vertical Tail")

## Assembly
wing_panels  = panel_wing(wing,                 # Wing or HalfWing type
                          [20, 3],              # Number of spanwise panels for half-wing 
                          10;                   # Chordwise panels 
                        #   spacing = "cosine"  # Spacing distribution: Default works well for symmetric
                         )

htail_panels = panel_wing(htail, [6], 6;
                          position	= [4., 0, 0],
                          angle 	= deg2rad(-2.),
                          axis 	  	= [0., 1., 0.],
                          spacing   = "uniform"
                         )
vtail_panels = panel_wing(vtail, [6], 5; 
                          position 	= [4., 0, 0],
                          angle 	= π/2, 
                          axis 	 	= [1., 0., 0.],
                          spacing   = "uniform"
                         )

aircraft = Dict("Wing" 			  	=> wing_panels,
                "Horizontal Tail" 	=> htail_panels,
                "Vertical Tail"   	=> vtail_panels)

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, α, β = 1.0, 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

@time data = 
    solve_case(aircraft, fs; 
               rho_ref     = ρ, 		# Reference density
               r_ref       = ref, 		# Reference point for moments
               area_ref    = S, 		# Reference area
               span_ref    = b, 		# Reference span
               chord_ref   = c, 		# Reference chord
               name        = ac_name,	# Aircraft name
               print       = true,		# Prints the results for the entire aircraft
            #    print_components = true,	# Prints the results for each component
              );

## Data collection
comp_names = (collect ∘ keys)(data) # Gets aircraft component names from analysis
comp  = comp_names[1]			   # Pick your component
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = data[comp]; #  Get the nearfield, farfield, force and moment coefficients, and other data for post-processing
print_coefficients(nf_coeffs, ff_coeffs, comp)

## Stability case
@time dv_data = 
    solve_stability_case(aircraft, fs;
                         rho_ref     = ρ,
                         r_ref       = ref,
                         area_ref    = S,
                         span_ref    = b,
                         chord_ref   = c,
                         name        = ac_name,
                         print       = true,
                         print_components = true,
                        );

## Data collection
names = (collect ∘ keys)(dv_data) 
comp  = names[1]
nf, ff, dvs = dv_data[comp];
print_coefficients(nf, ff, comp)
print_derivatives(dvs, comp)

## Streamlines

# Chordwise distribution
# using StaticArrays

# num_points = 100
# max_z = 2
# y = span(wing) / 10
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 10
init        = chop_leading_edge(wing, span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz]) ; 
                init .+ Ref([dx, dy,-dz]) ];

distance = 8
num_stream_points = 200
streams = plot_streams(fs, seed, horseshoes, Γs, distance, num_stream_points);

##
using Plots
gr(size = (200, 100), dpi = 300)

##
horseshoe_coords = plot_panels(horseshoe_panels[:]);

z_limit = b
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     camera = (30, 60),
    #  xlim = (-z_limit/2, z_limit/2),
     aspect_ratio = 1, 
     zlim = (-z_limit/2, z_limit/2),
     size = (1280, 720)
    )
plot!.(horseshoe_coords, color = :gray, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()

## Exaggerated CF distribution for plot, only works with GR and not Plotly

# Forces
wind_CFs = body_to_wind_axes.(CFs, fs.alpha, fs.beta)
CDis     = @. getindex(wind_CFs, 1)
CYs	     = @. getindex(wind_CFs, 2)
CLs      = @. getindex(wind_CFs, 3)

hs_pts = bound_leg_center.(horseshoes)
hs_xs  = getindex.(hs_pts, 1)
hs_ys  = getindex.(hs_pts, 2)
hs_zs  = getindex.(hs_pts, 3)

quiver!(hs_xs[:], hs_ys[:], hs_zs[:], quiver=(CDis[:], CYs[:], CLs[:]) .* 500)
plot!()