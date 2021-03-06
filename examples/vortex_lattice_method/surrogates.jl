## Example script for building surrogate models using VLM analysis, needs Surrogates and DataFrames installed
using Surrogates
using DataFrames
using AeroMDAO

## Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing       = Wing(foils     = wing_foils,
                  chords    = [1.0, 0.6],
                  twists    = [0.0, 0.0],
                  spans     = [5.0],
                  dihedrals = [11.3],
                  sweep_LEs = [2.29]);
wing_mac   = mean_aerodynamic_center(wing)
wing_pos   = [0., 0., 0.]
print_info(wing, "Wing")

# Horizontal tail
htail = Wing(chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39])
htail_mac = mean_aerodynamic_center(htail)
htail_pos = [5., 0., 0.]
α_h_i     = 0.
print_info(htail, "Horizontal Tail")

# Vertical tail
vtail_foil = Foil(naca4((0,0,0,9)))
vtail = HalfWing(foils     = fill(vtail_foil, 2), 
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97])
vtail_mac = mean_aerodynamic_center(vtail) # NEEDS FIXING FOR ROTATION
vtail_pos = [5., 0., 0.]
print_info(vtail, "Vertical Tail")

## Panelling and assembly
wing_panels  = panel_wing(wing, [20], 10,
                          position = wing_pos
                         );
htail_panels = panel_wing(htail, [12], 12;
                          position	= htail_pos,
                          angle 	= deg2rad(α_h_i),
                          axis 	  	= [0., 1., 0.]
                         )
vtail_panels = panel_wing(vtail, [12], 10; 
                          position 	= vtail_pos,
                          angle 	= π/2, 
                          axis 	 	= [1., 0., 0.]
                         )

aircraft     = Dict(
                    "Wing"            => wing_panels,
                    "Horizontal Tail" => htail_panels,
                    "Vertical Tail"   => vtail_panels
                   )

## VLM setup
function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
    # Evaluate case
    data =  solve_case(aircraft, fs; 
                       rho_ref   = ρ, 
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print     = print
                      );
    
    # Get data
    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = data["Aircraft"]

    # Filter relevant data
    [ ff_coeffs[1:3]; nf_coeffs[4:6] ]
end

## Evaluate one case
# Case
x_w     = wing_pos + [ wing_mac[1], 0, 0 ]
ac_name = "My Aircraft"
ρ       = 1.225
ref     = x_w
V, α, β = 1.0, 0.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

@time coeffs = vlm_analysis(aircraft, fs, ρ, ref, S, b, c, true)

## Looping
αs, βs   = float.(-4:4), float.(-3:3)
fses     = [ Freestream(V, α, β, Ω) for α in αs, β in βs ]

## Evaluate cases
results = vlm_analysis.(Ref(aircraft), fses, ρ, Ref(ref), S, b, c, true);

## Generate DataFrame
αβs  = [ (α, β) for α in αs, β in βs ]
data = DataFrame([ (α, β, CD, CY, CL, Cl, Cm, Cn) for ((α, β), (CD, CY, CL, Cl, Cm, Cn)) in zip(αβs, results) ][:])
rename!(data, [:α, :β, :CD, :CY, :CL, :Cl, :Cm, :Cn])

##
using Plots
plotly(size = (1280, 720), markersize = 2, dpi = 300)

## The labels don't work for some reason, need to fix that later...
CD_plot = scatter(data[:,1], data[:,2], data[:,3], label = "CD")
CY_plot = scatter(data[:,1], data[:,2], data[:,4], label = "CY")
CL_plot = scatter(data[:,1], data[:,2], data[:,5], label = "CL")
Cl_plot = scatter(data[:,1], data[:,2], data[:,6], label = "Cl")
Cm_plot = scatter(data[:,1], data[:,2], data[:,7], label = "Cm")
Cn_plot = scatter(data[:,1], data[:,2], data[:,8], label = "Cn")

plots = [ CD_plot CY_plot CL_plot Cl_plot Cm_plot Cn_plot ]
plot(plots..., layout = (2, 3), xlabel = "α, ᵒ", ylabel = "β, ᵒ", zlabel = :nothing)
plot!()

## Surrogate training data
lower_bound, upper_bound = minimum.(eachcol(data[:,3])), maximum.(eachcol(data[:,3]))
zs = data[:,3]

## Generate surrogate
surrogate_model = LobachevskySurrogate(αβs, zs, lower_bound, upper_bound); 

## Evaluate surrogate
num     = 100
xs, ys  = range(-4, 4, length = num), range(-3, 3, length = num)
xys     = [ (x, y) for x in xs, y in ys ]
surr_zs = surrogate_model.(xys)

## Create DataFrame
surr_data = DataFrame((x, y, z) for ((x, y), z) in zip(xys, surr_zs))
rename!(surr_data, [:α, :β, :CD])

## Plot surrogate
scatter(data[:,1], data[:,2], data[:,3], xlabel = "α, ᵒ", ylabel = "β, ᵒ", zlabel = "CD", label = "CD")
mesh3d!(surr_data[:,1], surr_data[:,2], surr_data[:,3], label = "Surrogate CD")