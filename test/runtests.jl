using AeroMDAO
using Test

@testset "Kulfan CST Doublet-Source Panel Method" begin
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

@testset "Airfoil Processing" begin
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

@testset "Trapezoidal Wing Geometry" begin
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

@testset "NACA-4 Wing Vortex Lattice Method" begin
    # Define wing
    wing = Wing(foils     = Foil.(naca4((2,4,1,2)) for i ∈ 1:3),
                chords    = [0.18, 0.16, 0.08],
                twists    = [2., 0., -2.],
                spans     = [0.5, 0.5],
                dihedrals = [0., 11.3],
                sweep_LEs = [1.14, 8.])

    # Define reference values
    ρ   = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω   = [0.0, 0.0, 0.0]
    
    # Define freestream condition
    uniform = Freestream(10.0, 2.0, 2.0, Ω)

    # Evaluate stability case
    nf_coeffs, ff_coeffs, dv_coeffs = solve_stability_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5)

    # Test values
    nf_tests = [0.00094437, -0.00049066, 0.11137544, -0.001868206, -0.00440614, -9.807e-5, 0.0, 0.0, 0.0]
    ff_tests = [0.00103027, -0.00053694, 0.11136510]
    dv_tests = [ 0.0232338    3.17931e-5    0.00233895  -0.0108045     0.000566369  ;
                 0.00307003  -0.0139215     0.157886     0.000211956  -0.0193726    ;
                 1.56531794   0.000698079  -0.0340142    0.414198      0.00533197   ;
                 0.00904046  -0.0529491     0.737322     0.00279386   -0.0872407    ;
                -0.0323204   -0.00384898    0.0267506   -0.0132855    -0.00322571   ;
                -0.00180434  -0.00281951    0.0306778   -0.000785897   0.00249871   ]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nf_coeffs, nf_tests) ]

    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ff_coeffs, ff_tests) ]

    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dv_coeffs, dv_tests) ]
end