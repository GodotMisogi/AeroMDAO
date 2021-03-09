## 
using Revise
using StaticArrays
using AeroMDAO
using BenchmarkTools

## Wing section setup
foil        =   naca4((2,4,1,2))
wing_right  =   HalfWing(Foil.(foil for i ∈ 1:3),   # Foils
						[0.18, 0.16, 0.08],         # Chords
						[2., 0., -2.],              # Twists
						[0.5, 0.2],                 # Spans
						[0., 11.3],                 # Dihedrals
						[1.14, 8.])                 # Sweeps

wing 		= Wing(wing_right, wing_right)

print_info(wing)

## Assembly
ρ 			= 1.225
ref 		= [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
V, α, β 	= 10.0, 5.0, 5.0
Ω 			= [0.0, 0.0, 0.0]
freestream 	= Freestream(V, α, β, Ω)
@time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, freestream, ρ, ref; span_num = 25, chord_num = 10) 

#
begin
	println("\nNearfield:")
	print_dynamics(nf_coeffs...)
	println("\nFarfield:")
	print_dynamics(ff_coeffs...)
end

## Plotting
using Plots

## Coordinates
horseshoe_coords 	= plot_panels(horseshoe_panels[:])
camber_coords		= plot_panels(camber_panels[:])
wing_coords 		= plot_wing(wing);

CDis = getindex.(CFs, 1)
CYs	 = getindex.(CFs, 2)
CLs  = getindex.(CFs, 3);

colpoints 	= horseshoe_point.(horseshoe_panels)
xs 			= getindex.(colpoints, 1);
ys 			= getindex.(colpoints, 2);
zs 			= getindex.(colpoints, 3);
cl_pts 		= tupvector(SVector.(xs[:], ys[:], zs[:] .+ CLs[:]));

## Streamlines

# num_points = 50
# max_z = 0.1
# y = span(wing) / 2 - 0.05
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

span_points = 20
init        = trailing_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz])  ; 
				init .+ Ref([dx, dy, -dz]) ];

distance = 2
num_stream_points = 100
streams = plot_streams(freestream, seed, horseshoes, Γs, distance, num_stream_points);

## Display
z_limit = 0.5
plot(xaxis = "x", yaxis = "y", zaxis = "z",
	 aspect_ratio = 1, 
	 camera = (30, 30),
	 zlim = (-0.1, z_limit),
	 size = (800, 600))
# plot!.(horseshoe_coords, color = :black, label = :none)
plot!.(camber_coords, color = :black, label = :none)
scatter!(tupvector(colpoints)[:], marker = 1, color = :black, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()

## Span forces
plot1 = plot(ys[1,:], sum(CDis, dims = 1)[:], label = :none, ylabel = "CDi")
plot2 = plot(ys[1,:], abs.(sum(CYs, dims = 1)[:]), label = :none, ylabel = "CY")
plot3 = plot(ys[1,:], sum(CLs, dims = 1)[:], label = :none, xlabel = "y", ylabel = "CL")
plot(plot1, plot2, plot3, layout = (3,1))

## Lift distribution
plot(xaxis = "x", yaxis = "y", zaxis = "z",
	aspect_ratio = 1,
	camera = (30, 30),
	zlim = (-0.1, z_limit))
plot!(wing_coords, label = :none)
scatter!(cl_pts, zcolor = CLs[:], marker = 1, label = "CL")
plot!(size = (800, 600))