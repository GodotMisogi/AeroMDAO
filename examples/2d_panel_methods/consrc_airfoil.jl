## Example script for constant-strength source panel method
using LinearAlgebra
using Base.Iterators
using Seaborn
using AeroMDAO

## Airfoil
airfoil = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = true))
V, α    = 1., 0.
ρ       = 1.225
uniform = Uniform2D(V, α)
num_pans = 80
panels  = paneller(airfoil, num_pans);

## Constant-strength source panel
A = @time constant_source_matrix(panels)
b = @time constant_source_boundary_condition(panels, velocity(uniform))
σs = A \ b

##
sum(σs)

## Pressure coefficient
pts = collocation_point.(panels);
panel_vels = [ velocity(uniform) .+ sum(source_velocity.(σs, panels, Ref(panel_i))) for panel_i in panels ] 

qts = @. dot(panel_vels, panel_tangent(panels))
cps = @. 1 - qts^2 / uniform.mag^2

##
upper, lower = get_surface_values(panels, cps)
lower = [ upper[end]; lower; upper[1] ]

plot(first.(upper), last.(upper), label = "Upper")
plot(first.(lower), last.(lower), label = "Lower")
ylim(maximum(cps), minimum(cps))
xlabel("(x/c)")
ylabel("Cp")
# legend()
show()

## Plotting
x_domain, y_domain = (-1, 2), (-1, 1)
grid_size = 50
x_dom, y_dom = range(x_domain..., length = grid_size), range(y_domain..., length = grid_size)
grid = product(x_dom, y_dom)
pts = panel_points(panels);

source_vels = [ velocity(uniform) .+ sum(source_velocity.(σs, panels, x, y)) for (x, y) in grid ]

speeds = @. sqrt(first(source_vels)^2 + last(source_vels)^2)
cps = @. 1 - speeds^2 / uniform.mag^2

contourf(first.(grid), last.(grid), cps)
CB1 = colorbar(label = "Pressure Coefficient (Cp)")
# quiver(first.(grid), last.(grid), first.(source_vels), last.(source_vels), speeds)
streamplot(first.(grid)', last.(grid)', first.(source_vels)', last.(source_vels)', color = speeds', cmap = "coolwarm", density = 3)
CB2 = colorbar(orientation="horizontal", label = "Relative Airspeed (U/U∞)")
fill(first.(pts), last.(pts), color = "black", zorder = 3)
tight_layout()
show()