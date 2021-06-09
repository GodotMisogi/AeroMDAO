using Revise
using Rotations
using AeroMDAO
using BenchmarkTools
using ForwardDiff
using StaticArrays

## Wing
TrapezoidalWing(b, δ, Λ, λ, c_root, τ_root, τ_tip) = HalfWing([c_root, λ * c_root],	[τ_root, τ_tip], [b], [δ], [Λ])

xs = [4.0, 0.0, 15.0, 0.4, 2.0, 0.0, -2.0]
wing_right  = TrapezoidalWing(xs...)
wing        = Wing(wing_right, wing_right)
wing_mac 	= mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Differentiation tests
winger(x) = let wing = TrapezoidalWing(x...); Wing(wing, wing) end
ForwardDiff.jacobian(mean_aerodynamic_center ∘ winger, xs)

## VLM
function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
	# Evaluate case
	data = 	solve_case(aircraft, fs; 
                       rho_ref   = ρ, 
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print 	 = print
                      );
		
	nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = data["Aircraft"]
			   
	[ ff_coeffs[1:3]; nf_coeffs[4:6] ]
end

function vlmer(x)
	wing 		= winger(x)
	wing_pos    = [0., 0., 0.]
	ρ 			= 1.225
	ref     	= zeros(3)
	V, α, β 	= 1.0, 0.0, 0.0
	Ω 			= zeros(3)
	fs 	    	= Freestream(V, α, β, Ω)
	S, b, c 	= projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
	wing_panels = panel_wing(wing, [20], 10,
                          	 position = wing_pos);
	aircraft 	= Dict("Wing" => wing_panels)
	results 	= vlm_analysis(aircraft, fs, ρ, ref, S, b, c)
end

vlmer(xs)

## ForwardDiff
ForwardDiff.jacobian(vlmer, [5.0, 0.0, 10.0, 0.4, 2.0, 0.0, -2.0])