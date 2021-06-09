## Mesh convergence tests
using AeroMDAO
using Plots
using DataFrames, StatsPlots

## Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 3))
wing_right = HalfWing(wing_foils,
                      [1.0, 0.6, 0.2],
                      [0.0, 0.0, 0.0],
                      [5.0, 0.5],
                      [5., 5.],
                      [5., 5.]);
wing = Wing(wing_right, wing_right)
print_info(wing, "Wing")
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Assembly
ρ       = 1.225
ref     = [0.25, 0., 0.]
V, α, β = 1.0, 5.0, 5.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

## Panelling error closure
panel_case(x, y) = solve_case(wing, fs; 
                              rho_ref   = ρ, 
                              r_ref     = ref,
                              area_ref  = S,
                              span_ref  = b,
                              chord_ref = c,
                              span_num  = y, 
                              chord_num = x,
                              viscous   = false,
                              x_tr      = [0.3, 0.3]);

## Distributions
xs = 6:30
ys = xs .* 4

nf_data, ff_data = zeros(5)', zeros(5)'
for (x, y) in zip(xs, ys)
    @time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = 
    panel_case(x, y)

    println("Distribution: nx = $x, ny = $y")	
    print_coefficients(nf_coeffs, ff_coeffs)
    nf_data = [ nf_data; [x; y; nf_coeffs[1:3] ]' ]
    ff_data = [ ff_data; [x; y; ff_coeffs[1:3] ]' ]
end

##
data = DataFrame([ nf_data ff_data[:,3:end] ], :auto)
rename!(data, [:chord_num, :span_num, :CD_nf, :CY_nf, :CL_nf, :CD_ff, :CY_ff, :CL_ff])

##
errors = @. abs((data[3:end,3:end] - data[2:end-1,3:end]) / data[2:end-1,3:end])

##
plot(xlabel = "Spanwise Panels", ylabel = "Error", yscale = :log10)
plot!(data[2:end-1,"span_num"], errors[!,"CD_nf"], label = "CD Nearfield")
plot!(data[2:end-1,"span_num"], errors[!,"CY_nf"], label = "CY Nearfield")
plot!(data[2:end-1,"span_num"], errors[!,"CL_nf"], label = "CL Nearfield")
plot!(data[2:end-1,"span_num"], errors[!,"CD_ff"], label = "CD Farfield")
plot!(data[2:end-1,"span_num"], errors[!,"CY_ff"], label = "CY Farfield")
plot!(data[2:end-1,"span_num"], errors[!,"CL_ff"], label = "CL Farfield")

##
plot()
plot!(data[2:end,"span_num"], data[2:end,"CL_nf"], marker = :dot, label = "Nearfield")
plot!(data[2:end,"span_num"], data[2:end,"CL_ff"], marker = :dot, label = "Farfield") 