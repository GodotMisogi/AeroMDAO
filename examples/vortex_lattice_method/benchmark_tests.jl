## Aircraft analysis benchmarking
using BenchmarkTools
using StaticArrays

# All subsequent analyses use no symmetry tricks for performance as AeroMDAO hasn't implemented it and apples must be compared to apples.

## AeroMDAO tests: https://github.com/GodotMisogi/AeroMDAO.jl
#=======================================================#
using AeroMDAO

println("AeroMDAO Aircraft Functional -")
@benchmark begin
    # Wing
    wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                chords    = [1.0, 0.6],
                twists    = [2.0, 2.0],
                spans     = [5.0],
                dihedrals = [11.31],
                sweep_LEs = [2.29]);

    # Horizontal tail
    htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.25],
                 dihedrals = [0.],
                 sweep_LEs = [6.39])

    # Vertical tail
    vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                     chords    = [0.7, 0.42],
                     twists    = [0.0, 0.0],
                     spans     = [1.0],
                     dihedrals = [0.],
                     sweep_LEs = [7.97])

    wing_panels  = panel_wing(wing, 20, 20, spacing = "cosine")

    htail_panels = panel_wing(htail, 12, 12;
                              position = [4., 0, 0],
                              angle    = deg2rad(0.),
                              axis 	   = [0., 1., 0.],
                              spacing  = "cosine"
                             )

    vtail_panels = panel_wing(vtail, 12, 10; 
                              position = [4., 0, 0],
                              angle    = π/2, 
                              axis 	   = [1., 0., 0.],
                              spacing  = "cosine"
                             )

    aircraft = Dict("Wing" 			  	=> wing_panels,
                    "Horizontal Tail" 	=> htail_panels,
                    "Vertical Tail"   	=> vtail_panels)

    # display(size.([ wing_panels[1], htail_panels[1], vtail_panels[1] ])) # Checking sizes

    ρ       = 1.225
    x_ref   = [0.5, 0., 0.]
    V, α, β = 1.0, 5.0, 0.0
    Ω       = [0.0, 0.0, 0.0]
    fs      = AeroMDAO.Freestream(V, α, β, Ω)
    S, b, c = 9.0, 10.0, 0.9

    data = solve_case(aircraft, fs; 
                      rho_ref   = ρ,
                      r_ref     = x_ref,
                      area_ref  = S,
                      span_ref  = b,
                      chord_ref = c)
    
    nf, ff = data["Aircraft"][1:2]

    nf[1:3], nf[4:6], ff[1]
end

## BYU FLOW Lab tests: https://github.com/byuflowlab/VortexLattice.jl
#=======================================================#

# Example taken from: https://flow.byu.edu/VortexLattice.jl/dev/examples/

# Pulled these out to avoid compilation for proper benchmarking
fc   = fill((xc) -> 0, 2) # camberline function for each section
fc_h = fill((xc) -> 0, 2) # camberline function for each section
fc_v = fill((xc) -> 0, 2) # camberline function for each section

##
using VortexLattice

println("BYU FLOW Lab VortexLattice.jl - ")
@benchmark begin
    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 20
    nc = 20
    spacing_s = Cosine()
    spacing_c = Cosine()
    mirror = true

    # horizontal stabilizer
    xle_h = [0.0, 0.14]
    yle_h = [0.0, 1.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.42]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    ns_h = 12
    nc_h = 12
    spacing_s_h = Cosine()
    spacing_c_h = Cosine()
    mirror_h = true

    # vertical stabilizer
    xle_v = [0.0, 0.14]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 1.0]
    chord_v = [0.7, 0.42]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    ns_v = 12
    nc_v = 10
    spacing_s_v = Cosine()
    spacing_c_v = Cosine()
    mirror_v = false

    Sref = 9.0
    cref = 0.9
    bref = 10.0
    rref = [0.5, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = VortexLattice.Freestream(Vinf, alpha, beta, Omega)

    symmetric = [false, false, false]

    # generate surface panels for wing
    wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [4.0, 0.0, 0.0])
    translate!(htail, [4.0, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, [4.0, 0.0, 0.0])
    translate!(vtail, [4.0, 0.0, 0.0])

    # display(size.([ wing, htail, vtail ])) # Checking sizes

    grids = [wgrid, hgrid, vgrid]
    surfaces = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=Wind())

    CDiff = far_field_drag(system)

    CF, CM, CDiff
end