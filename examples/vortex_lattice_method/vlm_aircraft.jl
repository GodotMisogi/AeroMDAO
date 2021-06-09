## Aircraft analysis case
using AeroMDAO

## Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing_right = HalfWing(wing_foils,
                      [1.0, 0.6],
                      [2.0, 2.0],
                      [5.0],
                      [11.3],
                      [2.29]);
wing = Wing(wing_right, wing_right)
print_info(wing, "Wing")

# Horizontal tail
htail_foil = Foil(naca4((0,0,1,2)))
htail_right = HalfWing(fill(htail_foil, 2),
                       [0.7, 0.42],
                       [0.0, 0.0],
                       [1.25],
                       [0.],
                       [6.39])
htail = Wing(htail_right, htail_right)
print_info(htail, "Horizontal Tail")

# Vertical tail
vtail_foil = Foil(naca4((0,0,0,9)))
vtail = HalfWing(fill(vtail_foil, 2), 
                      [0.7, 0.42],
                      [0.0, 0.0],
                      [1.0],
                      [0.],
                      [7.97])
print_info(vtail, "Vertical Tail")

# Assembly
wing_panels  = panel_wing(wing, [20], 10);
htail_panels = panel_wing(htail, [12], 12;
                          position	= [4., 0, 0],
                          angle 	= deg2rad(-2.),
                          axis 	  	= [0., 1., 0.]
                         )
vtail_panels = panel_wing(vtail, [12], 10; 
                          position 	= [4., 0, 0],
                          angle 	= π/2, 
                          axis 	 	= [1., 0., 0.]
                         )

aircraft = Dict("Wing" 			  	=> wing_panels,
                "Horizontal Tail" 	=> htail_panels,
                "Vertical Tail"   	=> vtail_panels)

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

# Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, α, β = 1.0, 1.0, 1.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

data = 
    solve_case(aircraft, fs; 
               rho_ref     = ρ, 		# Reference density
               r_ref       = ref, 		# Reference point for moments
               area_ref    = S, 		# Reference area
               span_ref    = b, 		# Reference span
               chord_ref   = c, 		# Reference chord
               name        = ac_name,	# Aircraft name
               print       = true,		# Prints the results for the entire aircraft
               print_components = true,	# Prints the results for each component
              );

## Data collection
names = (collect ∘ keys)(data) # Gets aircraft component names from analysis
comp  = names[1]			   # Pick your component
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = data[comp]; #  Get the nearfield, farfield, force and moment coefficients, and other data for post-processing
print_coefficients(comp, nf_coeffs, ff_coeffs)

## Stability case
dv_data = 
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
print_coefficients(comp, nf, ff)
print_derivatives(comp, dvs)

## Streamlines

# Chordwise distribution
# using StaticArrays

# num_points = 100
# max_z = 2
# y = span(wing) / 10
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 50
init        = leading_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz]) ; 
                init .+ Ref([dx, dy,-dz]) ];

distance = 8
num_stream_points = 200
streams = plot_streams(fs, seed, horseshoes, Γs, distance, num_stream_points);

## Plot
horseshoe_coords = plot_panels(horseshoe_panels[:])

##
using Plots
gr(size = (200, 100), dpi = 300)

##
z_limit = b
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1, 
     camera = (30, 60),
     xlim = (-z_limit/2, z_limit/2),
     zlim = (-z_limit/2, z_limit/2),
     size = (1280, 720)
    )
plot!.(horseshoe_coords, color = :black, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()