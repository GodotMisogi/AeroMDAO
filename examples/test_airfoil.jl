##
using Revise
using AeroMDAO

##
foilpath = joinpath((dirname ∘ dirname ∘ pathof)(AeroMDAO), "test/CRM.dat")

## Foil processing
coords = read_foil(foilpath)

## Cosine spacing
cos_foil = cosine_foil(coords, 51)

## Camber-thickness transformations
xcamthick = foil_camthick(cos_foil)
foiler = camthick_foil(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])

## Fitting to CST
num_dv = 4

# Coordinates fitting
up, low = split_foil(cos_foil)
alpha_u, alpha_l = coords_to_CST(up, num_dv), coords_to_CST(low, num_dv)
cst_foil = kulfan_CST(alpha_u, alpha_l, (1e-4, -1e-4), 0.0)

## Camber-thickness fitting
alphas = camthick_to_CST(cos_foil, num_dv)
cam_foil = camber_CST(alphas..., (0., 2e-4), 0)

## Kulfan CST
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
dzs = (1e-4, 1e-4)
foil = kulfan_CST(alpha_u, alpha_l, dzs, 0.2)

## NACA 4-digit airfoils
naca = naca4((2,4,1,2))

## Plotting library
using Plots
plotlyjs();

## Cosine and camber-thickness
plot(getindex.(cos_foil, 1), getindex.(cos_foil, 2), label = "Cosine")
plot!(xcamthick[:,1], xcamthick[:,2], label = "Camber")
plot!(xcamthick[:,1], xcamthick[:,3], label = "Thickness", aspectratio = 1)
plot!(getindex.(cst_foil, 1), getindex.(cst_foil, 2), label = "CST Coordinates Fit")
plot!(getindex.(cam_foil, 1), getindex.(cam_foil, 2), label = "CST Camber-Thickness Fit", aspectratio = 1)
plot!(getindex.(foil, 1), getindex.(foil, 2), label = "CST", aspectratio = 1)
plot!(getindex.(naca, 1), getindex.(naca, 2), label = "NACA", aspectratio = 1)