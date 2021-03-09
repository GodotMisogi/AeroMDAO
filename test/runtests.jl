using AeroMDAO
using Test

@testset "Kulfan CST Doublet-Source Panel Method" begin
    alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
    alpha_l = [-0.2, -0.1, -0.1, -0.001]
    dzs     = (1e-4, 1e-4)

    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
    
    uniform = Uniform2D(1., 5.)
    cl, cls, cps, panels = solve_case(airfoil, uniform, num_panels = 100)

    @test isapprox(cl, 0.90827913, atol = 1e-6)
    @test isapprox(sum(cls), 0.99739457, atol = 1e-6)
end

@testset "Airfoil Processing" begin
    foilpath = joinpath((dirname ∘ dirname ∘ pathof)(AeroMDAO), "test/CRM.dat")
    coords = read_foil(foilpath)

    cos_foil = cosine_foil(coords, 51)

    up, low = split_foil(cos_foil)

    num_dv = 4
    alpha_u, alpha_l = coords_to_CST(up, num_dv), coords_to_CST(low, num_dv)

    cst_foil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0.0)

    uniform = Uniform2D(1., 5.)
    cl, cls, cps, panels = solve_case(cst_foil, uniform, num_panels = 60)

    @test isapprox(cl, 0.97077114, atol = 1e-6)
    @test isapprox(sum(cls), 1.06741494, atol = 1e-6)
end

@testset "NACA-4 Vortex Lattice Method" begin
    foil        =   naca4((2,4,1,2))
    wing_right  =   HalfWing(Foil.(foil for i ∈ 1:3),
                            [0.18, 0.16, 0.08],
                            [2., 0., -2.],
                            [0.5, 0.5],
                            [0., 11.3],
                            [1.14, 8.])

    wing = Wing(wing_right, wing_right)

    ρ = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω = [0.0, 0.0, 0.0]
    uniform = Freestream(10.0, 5.0, 5.0, Ω)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ρ, ref, span_num = 20, chord_num = 5) 

    CL_nf, CDi_nf, CY_nf, Cl_nf, Cm_nf, Cn_nf, p_b_nf, q_b_nf, r_b_nf = nf_coeffs
    CL_ff, CDi_ff, CY_ff, Cl_ff, Cm_ff, Cn_ff, p_b_ff, q_b_ff, r_b_ff = ff_coeffs
    
    @test isapprox(CL_ff, 0.66478722, atol = 1e-6)
    @test isapprox(CDi_ff, 0.00099817, atol = 1e-6)
    @test isapprox(CY_ff, -0.00359138, atol = 1e-6)
    @test isapprox(Cl_nf, 0.00792433, atol = 1e-6)
    @test isapprox(Cm_nf, -0.13004060, atol = 1e-6)
    @test isapprox(Cn_nf, 0.00076269, atol = 1e-6)
end