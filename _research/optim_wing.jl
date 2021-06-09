## 
using Revise
using StaticArrays
using AeroMDAO
using TimerOutputs
using Optim
using PlotlyJS

# Wing setup
#============================================#

## Helper functions
function make_wing(digits, chords, twists, spans, dihedrals, sweeps)
    wing_right = HalfWing(Foil.(naca4(tuple(digits...)) for i ∈ 1:5), chords, twists, spans, dihedrals, sweeps)
    Wing(wing_right, wing_right)
end


function run_case(wing, V, α)
    uniform = Freestream(V, α, 0.0)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 10, chord_num = 10)

    ff_coeffs
end

# Objective function
optimize_cdi(digits, chords, twists, spans, dihedrals, sweeps, V, α) = run_case(make_wing(digits, chords, twists, spans, dihedrals, sweeps), V, α)[2]

# Lift constraint function
cons_lift(digits, chords, twists, spans, dihedrals, sweeps, V, α, ρ = 1.225) = 
    let wing = make_wing(digits, chords, twists, spans, dihedrals, sweeps); dynamic_pressure(ρ, V) * projected_area(wing) * run_case(wing, V, α)[1] end

## Test runs
#============================================#

# Parameters
digits = [2,4,1,2]
chords = [0.18, 0.16, 0.08, 0.04, 0.02]
twists = [0., 0., 0., 0., 0.]
spans = [0.2, 0.2, 0.2, 0.2]
dihedrals = [0., 0., 0., 0.]
sweeps = [0., 0., 0., 0.]
V = 10.
α = 0.

reset_timer!()

@timeit "Making Wing" wing = make_wing(digits, chords, twists, spans, dihedrals, sweeps)

@timeit "Run Case" ff_coeffs = @time run_case(wing, V, α)

print_timer()

println("\nFarfield:")
print_dynamics(ff_coeffs...)

## Optimisation
#============================================#

optimize_chords = TwiceDifferentiable(x -> optimize_cdi(digits, x, twists, spans, dihedrals, sweeps, V, α), chords, autodiff = :forward)

l_bound = fill(1e-12, length(chords))
u_bound = fill(Inf, length(chords))
bound_chords = TwiceDifferentiableConstraints(l_bound, u_bound)

lc = [1.]; uc = [1.]

function lifter!(c, x) 
    c = cons_lift(digits, x, twists, spans, dihedrals, sweeps, V, α)
end

lift_constraint = TwiceDifferentiableConstraints(lifter!, l_bound, u_bound, lc, uc, :forward)

reset_timer!()

@timeit "Optimizing" res_chord = optimize(optimize_chords,
                                          lift_constraint,
                                        #   bound_chords,
                                          chords,							# Initial value
                                          IPNewton(),
                                          # autodiff = :forward,
                                          Optim.Options(
                                                          # extended_trace = true,
                                                          show_trace = true
                                                          )
                                          )

##
print_timer()

##
layout = Layout(title = "Vortex Lattice",
                scene=attr(aspectratio=attr(x=1,y=1,z=1)),
                )

opt_wing = make_wing(digits, res_chord.minimizer, twists, spans, dihedrals, sweeps)
wing_trace = trace_surface(opt_wing)

plot([ 
		(trace for trace in wing_trace)...,
	],
	layout)