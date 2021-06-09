## Aircraft analysis case
using AeroMDAO
using BenchmarkTools
using ForwardDiff

## Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing_right = HalfWing(wing_foils,
                      [1.0, 0.6],
                      [2.0, 2.0],
                      [5.0],
                      [11.3],
                      [2.29]);
wing = Wing(wing_right, wing_right)

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
htail_right = HalfWing(htail_foils,
                       [0.7, 0.42],
                       [0.0, 0.0],
                       [1.25],
                       [0.],
                       [6.39])
htail = Wing(htail_right, htail_right)

# Vertical tail
vtail_foils = Foil.(fill(naca4((0,0,0,9)), 2))
vtail = HalfWing(vtail_foils, 
                 [0.7, 0.42],
                 [0.0, 0.0],
                 [1.0],
                 [0.],
                 [7.97]);

## Assembly
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

aircraft = Dict(
                "Wing" 			  	=> wing_panels,
                "Horizontal Tail" 	=> htail_panels,
                "Vertical Tail"   	=> vtail_panels
                );

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, α, β = 1.0, 1.0, 1.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω);

## Pure
println("Pure -")
@time begin
    data = 
    solve_case(aircraft, fs; 
               rho_ref     = ρ, 		# Reference density
               r_ref       = ref, 		# Reference point for moments
               area_ref    = S, 		# Reference area
               span_ref    = b, 		# Reference span
               chord_ref   = c, 		# Reference chord
               name        = ac_name,	# Aircraft name
            #    print       = true,		# Prints the results for the entire aircraft
               print_components = true,	# Prints the results for each component
              );
end;

## Data collection
names = (collect ∘ keys)(data) 
comp  = names[1]			   
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = data[comp];
print_coefficients(nf_coeffs, ff_coeffs, comp)

## Impure
println("Impure -")
@time begin
    state = VLMState(fs.V, fs.alpha, fs.beta, fs.omega, 
    rho_ref   = ρ,
    r_ref     = ref,
    area_ref  = S, 
    chord_ref = c, 
    span_ref  = b, 
    name      = ac_name);

    system, surfs, nf_t, ff_t = AeroMDAO.VortexLattice.solve_case!(aircraft, state);
end;

##
surf_names = (collect ∘ keys)(surfs) 
comp = surf_names[2]			   
nf_coeffs, ff_coeffs, surf = surfs[comp];
# print_coefficients(nf_coeffs, ff_coeffs, comp)
print_coefficients(nf_t, ff_t, state.name)

## Derivatives
function wtf(aircraft, fs)
    
    # Closure
    function get_derivatives(x)
        state = VLMState(fs.V, x[1], x[2], x[3:end], 
                         rho_ref   = ρ,
                         r_ref     = ref,
                         area_ref  = S, 
                         chord_ref = c, 
                         span_ref  = b, 
                         name      = ac_name);

        system, surfs, nf_t, ff_t = AeroMDAO.VortexLattice.solve_case!(aircraft, state);

        nf_t
    end

    jac = ForwardDiff.jacobian(get_derivatives, [fs.alpha, fs.beta, fs.omega... ])
end

dvs = wtf(aircraft, fs)

##
