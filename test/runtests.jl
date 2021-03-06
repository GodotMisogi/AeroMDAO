using AeroMDAO
using Test

@testset "2D Panel Method - Doublet-Source Kulfan CST" begin
    # Define Kulfan CST coefficients
    alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
    alpha_l = [-0.2, -0.1, -0.1, -0.001]
    dzs     = (0., 0.)

    # Define airfoil
    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);

    # Define uniform flow
    uniform = Uniform2D(1., 5.)

    # Evaluate case
    cl, cls, cms, cps, panels = solve_case(airfoil, uniform; num_panels = 80)

    @test cl       ≈  0.84188988 atol = 1e-6
    @test sum(cls) ≈  0.84073703 atol = 1e-6
    @test sum(cms) ≈ -0.26104277 atol = 1e-6
end

@testset "Geometry - Airfoil Processing" begin
    # Import and read airfoil coordinates
    foilpath = joinpath((dirname ∘ dirname ∘ pathof)(AeroMDAO), "test/CRM.dat")
    coords   = read_foil(foilpath)

    # Cosine spacing
    cos_foil = cosine_foil(coords, 51)

    # Split airfoil
    up, low  = split_foil(cos_foil)

    # Convert coordinates to Kulfan CST variables
    num_dv   = 4
    alpha_u, alpha_l = coords_to_CST(up, num_dv), coords_to_CST(low, num_dv)

    # Generate same airfoil using Kulfan CST parametrisation
    cst_foil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0.0)

    uniform  = Uniform2D(1., 5.)
    cl, cls, cms, cps, panels = solve_case(cst_foil, uniform; num_panels = 80)

    @test cl       ≈  0.85736965 atol = 1e-6
    @test sum(cls) ≈  0.85976886 atol = 1e-6
    @test sum(cms) ≈ -0.29766116 atol = 1e-6
end

@testset "Geometry - Two-Section Trapezoidal Wing" begin
    # Define wing
    wing_right = HalfWing(chords    = [1.0, 0.6, 0.2],
                          twists    = [2.0, 0.0, -0.2],
                          spans     = [5.0, 0.5],
                          dihedrals = [5., 5.],
                          sweep_LEs = [5., 5.]);

    # Get wing info
    b, S, c, AR = info(wing_right)
    λ           = taper_ratio(wing_right)
    wing_mac    = mean_aerodynamic_center(wing_right)

    @test b        ≈ 5.50000000                    atol = 1e-6
    @test S        ≈ 4.19939047                    atol = 1e-6
    @test c        ≈ 0.79841269                    atol = 1e-6
    @test AR       ≈ 7.20342634                    atol = 1e-6
    @test λ        ≈ 0.20000000                    atol = 1e-6
    @test wing_mac ≈ [0.42092866, 1.33432539, 0.0] atol = 1e-6
end

@testset "Vortex Lattice Method - NACA 0012 Rectangular Wing" begin
    # Define wing
    wing = Wing(foils     = Foil.(naca4((0,0,1,2)) for i ∈ 1:2),
                chords    = [0.18, 0.16],
                twists    = [0., 0.],
                spans     = [0.5,],
                dihedrals = [5.],
                sweep_LEs = [1.14])

    # Define reference values
    ρ   = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω   = [0.0, 0.0, 0.0]

    # Define freestream condition
    uniform = Freestream(10.0, 2.0, 2.0, Ω)

    # Evaluate stability case
    nf_coeffs, ff_coeffs, dv_coeffs = solve_stability_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5)

    # Test values
    nf_tests = [0.001189, -0.000228, 0.152203, -0.000242, -0.003486, -8.1e-5, 0.0, 0.0, 0.0]
    ff_tests = [0.00123,  -0.000271, 0.152198]
    dv_tests = [ 0.068444 -0.000046 -0.000711  0.023607  0.000337;
                 0.010867 -0.007536  0.129968  0.021929 -0.012086;
                 4.402229 -0.012973 -0.070654  6.833903  0.001999;
                 0.031877 -0.013083  0.460035  0.091216 -0.039146;
                -0.112285 -0.004631  0.105695 -0.852395 -0.007696;
                -0.002218 -0.002115  0.008263 -0.003817  0.001079]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nf_coeffs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ff_coeffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dv_coeffs, dv_tests) ]
end

@testset "Vortex Lattice Method - Vanilla Aircraft" begin
    ## Wing
    wing = Wing(foils = Foil.(fill(naca4((0,0,1,2)), 2)),
                chords    = [1.0, 0.6],
                twists    = [0.0, 0.0],
                spans     = [5.0],
                dihedrals = [11.39],
                sweep_LEs = [0.]);

    # Horizontal tail
    htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
    htail = Wing(foils     = htail_foils,
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.25],
                 dihedrals = [0.],
                 sweep_LEs = [6.39])

    # Vertical tail
    vtail_foils = Foil.(fill(naca4((0,0,0,9)), 2))
    vtail = HalfWing(foils     = vtail_foils,
                     chords    = [0.7, 0.42],
                     twists    = [0.0, 0.0],
                     spans     = [1.0],
                     dihedrals = [0.],
                     sweep_LEs = [7.97])

    ## Assembly
    wing_panels  = panel_wing(wing, 16, 10;
                              spacing = "cosine")
    htail_panels = panel_wing(htail, 6, 6;
                              position = [4., 0, 0],
                              angle    = deg2rad(-2.),
                              axis     = [0., 1., 0.],
                              spacing  = "cosine")
    vtail_panels = panel_wing(vtail, 5, 6;
                              position  = [4., 0, 0],
                              angle     = π/2,
                              axis      = [1., 0., 0.],
                              spacing   = "cosine")

    aircraft = Dict("Wing"            => wing_panels,
                    "Horizontal Tail" => htail_panels,
                    "Vertical Tail"   => vtail_panels)

    ## Reference quantities
    ac_name = "My Aircraft"
    S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
    ρ       = 1.225
    ref     = [0.25c, 0., 0.]
    V, α, β = 1.0, 1.0, 1.0
    Ω       = [0.0, 0.0, 0.0]
    fs      = Freestream(V, α, β, Ω)

    ## Stability case
    dv_data = solve_stability_case(aircraft, fs;
                                   rho_ref   = ρ,
                                   r_ref     = ref,
                                   area_ref  = S,
                                   span_ref  = b,
                                   chord_ref = c,
                                   name      = ac_name);

    nfs, ffs, dvs = dv_data[ac_name]

    nf_tests = [0.000258, -0.006642, 0.074301, -0.003435, 0.075511, 0.001563, 0.0, 0.0, 0.0]
    ff_tests = [0.000375, -0.006685, 0.074281]
    dv_tests = [ 0.016795  0.003460  0.003761   0.093303 -0.000674;
                -0.000863 -0.374410  0.403476   0.000630 -0.253848;
                 5.749765  0.046649 -0.01346   15.571205  0.020396;
                 0.022674 -0.196605  0.660392   0.099065 -0.039688;
                -2.70367  -0.132928  0.070111 -37.372278 -0.064439;
                 0.002034  0.087382  0.014991   0.005840  0.091088]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end;
