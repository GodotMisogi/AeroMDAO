##
using AeroMDAO
using NLsolve

## Stability analysis
function aircraft_nearfield_forces(aircraft, fs, ρ, ref, S, b, c)
    nf, ff, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = 
        solve_case(aircraft, fs; 
                   rho_ref     = ρ, 
                   r_ref       = ref, 
                   area_ref    = S, 
                   span_ref    = b, 
                   chord_ref   = c,
                   )["Aircraft"]
    
    nf
end

pitching_moment_coefficient(aircraft, α :: Real, ρ, ref, S, b, c) = aircraft_nearfield_forces(aircraft, Freestream(1.0, α, 0., zeros(3)), ρ, ref, S, b, c)[5]

## Define aircraft
TrapezoidalWing(b, δ, Λ, λ, c_root, τ_root, τ_tip, foil_root, foil_tip) =
    HalfWing([ Foil(foil_root), Foil(foil_tip) ], # Foils
               [c_root, λ * c_root], 			  # Chords
               [τ_root, τ_tip], 				  # Twists
               [b],             				  # Span
               [δ],             				  # Dihedral
               [Λ])             				  # LE sweep

# Wing
wing_right  = TrapezoidalWing(4.0, 0.0, 15.0, 0.4, 2.0, 0.0, -2.0, naca4((2,4,1,2)), naca4((2,4,1,2)))
wing        = Wing(wing_right, wing_right)
wing_pos    = [0., 0., 0.];
S, b, c     = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
x_w         = wing_pos + [ c, 0, 0 ]

# Horizontal Tail
htail_right = TrapezoidalWing(1.0, 0.0, 15.0, 0.6, 0.8, 0.0, 0.0, naca4((0,0,1,2)), naca4((0,0,0,9)));
htail		= Wing(htail_right, htail_right)
htail_pos	= [5., 0., 0.]
α_h_i		= 0.;

# Vertical Tail
vtail		= TrapezoidalWing(0.8, 0.0, 8.0, 0.6, 0.8, 0.0, 0., naca4((0,0,0,9)), naca4((0,0,0,9)))
vtail_pos	= [5., 0., 0.];
                        
## Panelling and assembly
wing_panels  = panel_wing(wing, [20], 10;
                          position = wing_pos
                         )
htail_panels = panel_wing(htail, [10], 5;
                          position = htail_pos,
                          angle    = deg2rad(α_h_i),
                          axis 	   = [0., 1., 0.]
                         )
vtail_panels = panel_wing(vtail, [10], 5;
                          position = vtail_pos,
                          angle    = π/2
                         )

aircraft     = Dict("Wing" 			  => wing_panels,
                    "Horizontal Tail" => htail_panels,
                    "Vertical Tail"   => vtail_panels);

## Evaluate case
ac_name = "My Aircraft"
ρ 		= 1.225
ref     = x_w
V, α, β = 1.0, 0.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 	    = Freestream(V, α, β, Ω)

## Evaluate trim
trim_state!(R, α) = R -= pitching_moment_coefficient(aircraft, α[1], ρ, ref, S, b, c)
prob = nlsolve(trim_state!, [-1.0], 
               method = :newton,
            #    autodiff = :forward,
            #    show_trace = true, 
            #    extended_trace = true
               )
