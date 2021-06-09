## 
using Revise
using JuMP
using Ipopt
using StaticArrays
using AeroMDAO
using PlotlyJS
using ForwardDiff


# Wing setup
#============================================#

## Helper functions
function make_wing(x)
    wing_right = HalfWing(Foil.(naca4((0,0,1,2)) for i ∈ 1:5), 
                         [ x[1:5]... ], 
                         [0., 0., 0., 0., 0.],
                         [0.2, 0.2, 0.2, 0.2],
                         [0., 0., 0., 0.],
                         [ x[6:end]... ])
    span(Wing(wing_right, wing_right))
end

wing_chords = [0.18, 0.16, 0.08, 0.04, 0.02]
wing_sweeps = [20, 15, 10, 5.]
x0 = [ wing_chords; wing_sweeps ]

ForwardDiff.gradient(make_wing, x0)

##
function run_case(wing, V, α, β)
    ρ = 1.225
    Ω = SVector(0.0, 0.0, 0.0)
    fs = Freestream(V, α, β, Ω)
    c = mean_aerodynamic_chord(wing)
    ref = SVector(0.25 * c, 0., 0.)
    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, fs; rho_ref = ρ, r_ref = ref, area_ref = projected_area(wing), span_ref = span(wing), chord_ref = c, span_num = 80, chord_num = 20);
end

# Objective function
function optimize_cdi_chords(x...)
    wing = make_wing(x)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = run_case(wing, 10., 0., 0.)

    cdi = ff_coeffs[1]
end

function optimize_cdi_naca4(m, p, t, c)
    wing = naca4_wing(m, p, t, c)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = run_case(wing, 10., 0., 0.)
    
    cdi = ff_coeffs[1]
end

# Lift constraint function
function chords_lift(x...)
    ρ, speed = 1.225, 10.
    wing = make_wing(x)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = run_case(wing, speed, 0., 0.)
    
    cl = ff_coeffs[3]
    lift = 1/2 * ρ * speed^2 * projected_area(wing) * cl
end

function naca4_lift(m, p, t, c)
    ρ, speed = 1.225, 10.
    wing = naca4_wing(m, p, t, c)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = run_case(wing, speed, 0., 0.)

    cl = ff_coeffs[3]
    lift = 1/2 * ρ * speed^2 * projected_area(wing) * cl
end

## Test runs
#============================================#

wing_chords = [0.18, 0.16, 0.08, 0.04, 0.02]
wing_sweeps = [20, 15, 10, 5.]
x0 = [ wing_chords; wing_sweeps ]
wing = make_wing(x0)
cdi = optimize_cdi_chords(x0...)
lifter = chords_lift(x0...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

##
naca_test = (2,4,1,2)
wing = naca4_wing(naca_test...)
cdi = optimize_cdi_naca4(naca_test...)
lifter = naca4_lift(naca_test...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

## Chord optimization
#============================================#

chord_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 10

@variable(chord_design, 0.1 <= chordsweeps[1:num_dv - 1])

##
register(chord_design, :optimize_cdi_chords, num_dv - 1, optimize_cdi_chords, autodiff = true)

register(chord_design, :chords_lift, num_dv - 1, chords_lift, autodiff = true)

##
@NLobjective(chord_design, Min, optimize_cdi_chords(chordsweeps...))
@NLconstraint(chord_design, chords_lift(chordsweeps...) == 1)

## Run optimization
optimize!(chord_design)

## Print optimal case
println("Chords: $(value.(chords)), Optimal CDi: $(objective_value(chord_design))")
chord_vals = value.(chords)
uniform = Freestream(10., 0., 0.)
wing = make_wing(chord_vals...)
nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = run_case(wing, 10., 0., 0.)
cdi = ff_coeffs[2]
cl = ff_coeffs[1]
begin
    println("\nNearfield:")
    print_dynamics(nf_coeffs...)
    println("\nFarfield:")
    print_dynamics(ff_coeffs...)
end

##
layout = Layout(title = "Vortex Lattice",
                scene = attr(aspectratio=attr(x=1,y=1,z=1)));

## Streamlines
num_points = 30
max_z = 0.02
y = span(wing) / 2 - 0.5
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

span_points = 10
init = trailing_chopper(wing.right, span_points) 
dx, dy, dz = 0, 0, 1e-3
seed = [ init .+ Ref([dx, dy, dz]); 
         init .+ Ref([dx, dy, -dz]) ]

trace_streams = trace_streamlines(uniform, seed, horseshoes[:], Γs[:], 2, 100);

##
trace_horsies = trace_panels(horseshoe_panels[:])
trace_horses = trace_panels(horseshoe_panels[:], Γs[:])
trace_cambers = trace_panels(camber_panels[:])
trace_wing = trace_surface(wing)

PlotlyJS.plot(
            [
                (trace for trace in trace_horsies)...,
                (trace for trace in trace_horses)...,
                (trace for trace in trace_streams)...,
                (trace for trace in trace_cambers)...,
                # [ trace for trace in trace_wing ]...,
            ],
            layout)

## NACA-4 optimization
#============================================#

naca4_design = Model(with_optimizer(Ipopt.Optimizer))

@variable(naca4_design, 1e-4 <= digits[1:4] <= 1.) 

register(naca4_design, :optimize_cdi_naca4, 4, optimize_cdi_naca4, autodiff = true)

register(naca4_design, :naca4_lift, 4, naca4_lift, autodiff = true)

@NLobjective(naca4_design, Min, optimize_cdi_naca4(digits...))
@NLconstraint(naca4_design, naca4_lift(digits...) == 5)

## Run optimization
optimize!(naca4_design)

## Print
println("Digits: $(value.(chords)), Optimal CDi: $(objective_value(naca4_design))")
digit_vals = value.(digits)
cdi = optimize_cdi_chords(digit_vals...)
lifter = lift(digit_vals...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")