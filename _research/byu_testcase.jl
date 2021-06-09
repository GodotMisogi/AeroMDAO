using Test
using VortexLattice
using LinearAlgebra
using BenchmarkTools
using TimerOutputs

ztol = sqrt(eps())

# reflects a vector across the x-z plane
flipy(x) = [x[1], -x[2], x[3]]

# constructs a normal vector the way AVL does
function avl_normal_vector(dr, theta)
    st, ct = sincos(theta)
    bhat = dr/norm(dr) # bound vortex vector
    shat = [0, -dr[3], dr[2]]/sqrt(dr[2]^2+dr[3]^2) # chordwise strip normal vector
    chat = [ct, -st*shat[2], -st*shat[3]] # camberline vector
    ncp = cross(chat, dr) # normal vector perpindicular to camberline and bound vortex
    return ncp / norm(ncp) # normal vector used by AVL
end



# Wing and Tail without Finite Core Model

# NOTE: AVL's finite-core model is turned off for these tests

# NOTE: There is some interaction between twist, dihedral, and chordwise
# position which causes the normal vectors found by AVL to differ from those
# computed by this package.  We therefore manually overwrite the normal
# vectors when this occurs in order to get a better comparison.

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
ns = 40
nc = 20
spacing_s = Cosine()
spacing_c = Cosine()
mirror = false

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
mirror_h = false

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

# adjust chord lengths to match AVL (which uses chord length in the x-direction)
chord = @. chord/cos(theta)
chord_h = @. chord_h/cos(theta_h)
chord_v = @. chord_v/cos(theta_v)

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

symmetric = [true, true, false]

# reset_timer!()
@time begin
    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings - finite core deactivated
    wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        # @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = set_normal(p, ncp)
    end

    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    surfaces = [wing, htail, vtail]
    surface_id = [1, 1, 1]

    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    CDiff, CD, CY, CL, Cl, Cm, Cn
end
    # @test isapprox(CL, 0.60408, atol=1e-3)
    # @test isapprox(CD, 0.01058, atol=1e-4)
    # @test isapprox(CDiff, 0.010378, atol=1e-4)
    # @test isapprox(Cm, -0.02778, atol=2e-3)
    # @test isapprox(CY, 0.0, atol=ztol)
    # @test isapprox(Cl, 0.0, atol=ztol)
    # @test isapprox(Cn, 0.0, atol=ztol)

# print_timer()