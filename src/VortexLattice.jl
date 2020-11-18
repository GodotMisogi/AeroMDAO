module VortexLattice

include("MathTools.jl")

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs
using .MathTools: accumap, structtolist

export Laplace, Uniform3D, Uniform, velocity, 
Panel, Panel3D, area, panel_coords, midpoint, panel_normal,
stability_axes, wind_axes,
solve_horseshoes, dynamic_computations, nearfield_drag, 
streamlines, aerodynamic_coefficients, print_dynamics, horseshoe_lines

"""
Solutions to Laplace's equation.
"""
abstract type Laplace end

"""
A Uniform3D type expressing a uniform flow in spherical polar coordinates with angles α, β and magnitude mag.
"""
struct Uniform3D <: Laplace
    mag :: Real
    α :: Real 
    β :: Real
end

"""
Uniform3D constructor in degrees.
"""
Uniform(mag, α_deg, β_deg) = Uniform3D(-mag, deg2rad(α_deg), deg2rad(β_deg))

"""
Converts uniform flow coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r .* SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
Computes the velocity of Uniform3D.
"""
velocity(uni :: Uniform3D) = freestream_to_cartesian(uni.mag, uni.α, uni.β)

# Axis transformations
#==========================================================================================#

"""
Converts coordinates into stability axes.
"""
stability_axes(coords, uni :: Uniform3D) =  [cos(uni.α) 0 sin(uni.α); 
                                                  0     1     0     ;
                                             sin(uni.α) 0 sin(uni.α)] * coords
                                            # RotY{Float64}(uni.α) * coords

"""
Converts coordinates into wind axes. This is supposedly not used and is just ceremonially implemented.
"""
wind_axes(coords, uni :: Uniform3D) = RotZY{Float64}(uni.β, uni.α) * coords


"""
Flips the y axis of a given vector.
"""
y_flip(vector :: SVector{3, Float64}) = SVector(vector[1], -vector[2], vector[3])

"""
Flips the x and z axes of a given vector.
"""
stab_flip(vector :: SVector{3, Float64}) = SVector(-vector[1], vector[2], -vector[3])

"""
Transforms forces and moments into stability axes.
"""
stability_axes(force, moment, Ω, uniform :: Uniform3D) = stability_axes(force, uniform), stability_axes(stab_flip(moment), uniform), stability_axes(stab_flip(Ω), uniform)

"""
Transforms forces and moments into wind axes.
"""
wind_axes(force, moment, uniform :: Uniform3D) = wind_axes(force, uniform), wind_axes([ -moment[1], moment[2], -moment[3] ], uniform)


# Panel setup
#==========================================================================================#

"""
Placeholder. Panels should be an abstract type as they have some common methods, at least in 2D and 3D. 
"""
abstract type Panel end

"""
A composite type consisting of 4 coordinates. The following ASCII art depicts the order:

z → y
↓
x
        p1 —→— p4
        |       |
        ↓       ↓
        |       |
        p2 —→— p3
"""
struct Panel3D <: Panel
    p1 :: SVector{3,Real}
    p2 :: SVector{3,Real}
    p3 :: SVector{3,Real}
    p4 :: SVector{3,Real}
end

"""
Computes the area of Panel3D.
"""
area(panel :: Panel3D) = (abs ∘ norm)((panel.p2 - panel.p1) .* (panel.p3 - panel.p2))

"""
Computes the coordinates of a Panel3D for plotting purposes.
"""
panel_coords(panel :: Panel3D) = structtolist(panel)


"""
Computes the midpoint of Panel3D.
"""
midpoint(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) / 4

"""
Computes the normal vector of Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = panel.p3 .- panel.p1, 
                                     p42 = panel.p4 .- panel.p2, 
                                     p31_x_p42 = p31 × p42;
                                     p31_x_p42 / norm(p31_x_p42) end

"""
Computes the weighted value between two points.
"""
weighted(x1, x2, wx) = (1 - wx) * x1 + wx * x2

"""
Computes the weighted point between two points in a given direction.
"""
weighted_point((x1, y1, z1), (x2, y2, z2), wx, wy, wz) = SVector(weighted(x1, x2, wx), weighted(y1, y2, wy), weighted(z1, z2, wz))

"""
Computes the quarter point between two points in the x-z plane.
"""
quarter_point(p1, p2) = weighted_point(p1, p2, 1/4, 0, 1/4)

"""
Computes the 3-quarter point between two points in the x-z plane.
"""
three_quarter_point(p1, p2) = weighted_point(p1, p2, 3/4, 0, 3/4)

"""
Helper function to compute the bound leg of Panel3D for horseshoes/vortex rings, which is the quarter point on each side in the x-z plane.
"""
bound_leg(p1, p2, p3, p4) = [ quarter_point(p1, p2), quarter_point(p4, p3) ]

"""
Helper function to compute the collocation point of Panel3D for horseshoes/vortex rings, which is the 3-quarter point on each side in the x-z plane.
"""
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) .+ three_quarter_point(p4, p3) ) ./ 2

"""
A composite type consisting of two vectors to define a line.
"""
struct Line
    r1 :: SVector{3,Real}
    r2 :: SVector{3,Real}
end

"""
Computes the velocity induced at a point `r` by a vortex Line with uniform strength Γ. Checks if `r` lies on the line and sets it to `(0, 0, 0)` if so, as the velocity is singular there.
"""
function velocity(r, line :: Line, Γ, ε = 1e-8; trailing = false, direction = [1, 0, 0])
    a = r .- line.r1
    b = r .- line.r2
    a_x_b = a × b
    na, nb = norm(a), norm(b)

    # Compute bound leg velocity
    @timeit "Bound Leg" bound_velocity = any(<(ε), norm.([a, b, a_x_b])) ? zeros(3) : (1/na + 1/nb) * a_x_b / (na * nb + dot(a, b))
    
    # Compute velocities of trailing legs
    @timeit "Trailing Leg" trailing_velocity = !trailing ? zeros(3) : (a × direction) / (na - dot(a, direction)) / na - (b × direction) / (nb - dot(b, direction)) / nb

    # Sums bound and trailing legs' velocities
    @timeit "Sum Legs" Γ/(4π) * (bound_velocity .+ trailing_velocity)
end


"""
Placeholder. Vortex rings and horseshoes basically have the same methods, and are arrays of vortex lines.
"""
abstract type AbstractVortexArray end

"""
A horseshoe type consisting of a bound leg of type Line represening a vortex line.
"""
struct Horseshoe <: AbstractVortexArray
    bound_leg :: Line
end

"""
Computes the bound leg for a Panel3D.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Computes the collocation point of a Panel3D.
"""
collocation_point(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Returns a Horseshoe on a Panel3D.
"""
horseshoe_lines(panel :: Panel3D) = (Horseshoe ∘ Line)(bound_leg(panel)...)

"""
A vortex ring type consisting of vortex lines. TODO: Consider if better fields are "bound_vortices" and "trailing_vortices"
"""
struct VortexRing <: AbstractVortexArray
    left_leg :: Line
    bound_leg :: Line
    back_leg :: Line
    right_leg :: Line
end

"""
Helper function to compute the vortex ring given four points following Panel3D ordering.
"""
function vortex_ring(p1, p2, p3, p4)
    v1, v4 = quarter_point(p1, p2), quarter_point(p4, p3)
    v2, v3 = v1 .+ p2 .- p1, v4 .+ p3 .- p4
    [ v1, v2, v3, v4 ]
end

"""
Computes the vortex rings on a Panel3D.
"""
vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Constructor for vortex rings on a Panel3D using Lines. The following convention is adopted:

``` 
    p1 —bound_leg→ p4
    |               |
left_line       right_line
    ↓               ↓
    p2 —back_line→ p3
```
"""
function VortexRing(panel :: Panel3D)
    v1, v2, v3, v4 = vortex_ring(panel)

    bound_leg = Line(v1, v4)
    left_leg = Line(v2, v1)
    right_leg = Line(v3, v2)
    back_leg = Line(v4, v3)

    VortexRing(left_leg, bound_leg, back_leg, right_leg)
end

"""
Computes the midpoint of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_center(vortex :: AbstractVortexArray) = (vortex.bound_leg.r1 .+ vortex.bound_leg.r2) / 2

"""
Computes the direction vector of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_vector(vortex :: AbstractVortexArray) = vortex.bound_leg.r2 .- vortex.bound_leg.r1

"""
Computes the induced velocities at a point `r` of a Horseshoe with uniform strength Γ and trailing legs pointing in a given direction.
"""
velocity(r, horseshoe :: Horseshoe, Γ, norm_vel) = velocity(r, horseshoe.bound_leg, Γ, trailing = true, direction = norm_vel)

"""
Sums the velocities evaluated at a point `r` of vortex lines with uniform strength Γ.
"""
sum_vortices(r, vortex_lines :: Array{Line}, Γ) = sum(velocity(r, line, Γ) for line ∈ vortex_lines)

"""
Computes the induced velocities at a point `r` of a Vortex Ring with uniform strength Γ.
"""
velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)

# Matrix setup
#==========================================================================================#

"""
Computes the influence coefficient of the velocity of a vortex line at a collocation point projected to a normal vector.
"""
function influence_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal, norm_vel, symmetry = false) 
    if symmetry
        mirror_point = y_flip(collocation_point)
        col_vel = velocity(collocation_point, horseshoe, 1., norm_vel)
        mir_vel = (y_flip ∘ velocity)(mirror_point, horseshoe, 1., norm_vel)
        sum_vel = col_vel .+ mir_vel
        return dot(sum_vel, panel_normal)
    else 
        return dot(velocity(collocation_point, horseshoe, 1., norm_vel), panel_normal)
    end
end

"""
Assembles the Aerodynamic Influence Coefficient (AIC) matrix given vortex lines, collocation points and associated normal vectors.
"""
influence_matrix(colpoints_normals, vortex_lines :: Array{Horseshoe}, norm_vel, symmetry = false) = [ influence_coefficient(colpoint_i, horsie_j, normal_i, norm_vel, symmetry) for (colpoint_i, normal_i) ∈ colpoints_normals, horsie_j ∈ vortex_lines ]

"""
Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = @timeit "Dotting" dot.(velocities, normals)

"""
Generates a wake Panel3D aligned with a normalised direction vector.
"""
function wake_panel(panel :: Panel3D, uniform :: Uniform3D, bound = 1.)
    wake = bound .* velocity(uniform) / uniform.mag
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

"""
Generates wake Panel3Ds given the Panel3Ds of a wing corresponding to the trailing edge and a Uniform3D. The wake length and number of wake panels must be provided. 
"""
make_wake(last_panels :: Array{Panel3D}, uniform :: Uniform3D, wake_length, wake_num) = (permutedims ∘ accumap)(panel -> wake_panel(panel, uniform, wake_length / wake_num), wake_num, last_panels)

"""
Solves the AIC matrix with the boundary condition given Panel3Ds and a freestream velocity.
"""
function solve_horseshoes(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, norm_vel, Ω̄, symmetry = false) 

    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = [ norm_vel .+ (Ω̄ × rc_i) for rc_i ∈ colpoints ]

    
    @timeit "AIC" AIC = influence_matrix(zip(colpoints[:], normals[:]), horseshoes[:], norm_vel, symmetry) 
    @timeit "RHS" boco = boundary_condition(normals[:], total_vel[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = reshape(Γs, (length(horseshoes[:,1]), length(horseshoes[1,:]))), horseshoes

    output
end

# Force evaluations
#==========================================================================================#

# Trefftz plane evaluations
# downwash(xi, xj, zi, zj) = -1/(2π) * (xj - xi) / ( (zj - zi)^2 + (xj - xi)^2 )
# trefftz_plane(horseshoes :: Array{Horseshoes}) = [ [ horseshoe.vortex_lines[1].r1 for horseshoe ∈ horseshoes ]; 
#                                                    [ horseshoe.vortex_lines[3].r2 for horseshoe ∈ horseshoes ] ]
# function downwash(horseshoes :: Array{Horseshoes}) 
#     coords = trefftz_plane(horseshoes)
#     carts = product(coords, coords)
#     [ [ i == j ? 0 : downwash(x1, x2, z1, z2) for (i, (x1,y1,z1)) ∈ enumerate(coords) ] for (j, (x2,y2,z2)) ∈ enumerate(coords) ] 
# end

"""
Evaluates the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the aircraft reference frame.
"""
velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, freestream) = sum(velocity(r, horseshoe, Γ, freestream / norm(freestream)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .- freestream .- Ω × r


"""
Evaluates the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the global reference frame.
"""
velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, freestream) = sum(velocity(r, horseshoe, Γ, freestream / norm(freestream)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .+ freestream .+ Ω × r

"""
Computes the local Kutta-Jowkowski forces given an array of horseshoes, their associated vortex strengths Γs, a Uniform3D, and a density ρ. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
function kutta_force(Γs, horseshoes :: Array{Horseshoe}, vel, Ω, ρ)
    Γ_rings = zip(Γs, horseshoes)
    [ 
      let r_i = bound_leg_center(horseshoe_i), l_i = bound_leg_vector(horseshoe_i); 
      ρ * velocity(r_i, Ω, horseshoes, Γs, vel) × l_i * Γ_focus; 
      end for (Γ_focus, horseshoe_i) ∈ Γ_rings 
    ]
end

"""
Computes the dynamic pressure given density ρ and speed V.
"""
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: Array{Horseshoe}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) ∈ zip(forces, horseshoes) ]

"""
Computes the non-dimensional force coefficient corresponding to standard aerodynamics.
"""
force_coefficient(force, q, S) = force / (q * S)

"""
Computes the non-dimensional moment coefficient corresponding to standard flight dynamics.
"""
moment_coefficient(moment, q, S, ref_length) = force_coefficient(moment, q, S) / ref_length

"""
Computes the non-dimensional angular velocity coefficient corresponding to standard flight dynamics.
"""
rate_coefficient(angular_speed, V, ref_length) = angular_speed * ref_length / 2V

"""
Computes the Kutta-Jowkowski forces and associated moments for horseshoes.
"""
function dynamic_computations(Γs :: Array{<: Real}, horseshoes :: Array{<: AbstractVortexArray}, vel, Ω, r_ref, ρ = 1.225)
    @timeit "Forces" geom_forces = kutta_force(Γs, horseshoes, vel, Ω, ρ)
    @timeit "Moments" geom_moments = moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end

"""
Computes the near-field drag given the sum of the local Kutta-Jowkowski forces computed and a Uniform3D.
"""
nearfield_drag(force, uniform :: Uniform3D) = - dot(force, velocity(uniform) / uniform.mag)


# Post-processing
#==========================================================================================#

"""
Computes the streamlines from a given starting point, a Uniform3D, Horseshoes and their associated strengths Γs. The length of the streamline and the number of evaluation points must also be specified.
"""
function streamlines(point, uniform :: Uniform3D, Ω, horseshoes, Γs, length, num_steps)
    streamlines = [point]
    vel = velocity(uniform)
    for i ∈ 1:num_steps
        update = velocity(streamlines[end], Ω, horseshoes, Γs, vel)
        streamlines = [ streamlines..., streamlines[end] .+ (update / norm(update) * length / num_steps) ]
    end
    streamlines
end

"""
Computes the streamlines from the collocation points of given Horseshoes with the relevant previous inputs.
"""
streamlines(uniform :: Uniform3D, Ω, horseshoe_panels :: Array{<: Panel}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, length :: Real, num_steps :: Integer) = [ streamlines(SVector(hs), uniform, Ω, horseshoes, Γs, length, num_steps) for hs ∈ collocation_point.(horseshoe_panels)[:] ]

"""
Prints the relevant aerodynamic/flight dynamics information.
"""
function aerodynamic_coefficients(force, moment, drag, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)
    CDi = force_coefficient(drag, q, S)
    CY = force_coefficient(force[2], q, S)
    CL = force_coefficient(force[3], q, S)

    Cl = moment_coefficient(moment[1], q, S, b)
    Cm = moment_coefficient(moment[2], q, S, c)
    Cn = moment_coefficient(moment[3], q, S, b)

    p̄ = rate_coefficient(Ω[1], V, b)
    q̄ = rate_coefficient(Ω[2], V, c)
    r̄ = rate_coefficient(Ω[3], V, b)

    CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄
end

function print_dynamics(CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄)
    # println("Lift Coefficient")
    println("CL: $CL")
    # println("Drag Coefficient")
    println("CDi: $CDi")
    # println("Side Force Coefficient")
    println("CY: $CY")
    # println("Lift-to-Drag Ratio")
    println("L/D: $(CL/CDi)")
    # println("Rolling Moment Coefficient")
    println("Cl: $Cl")
    # println("Pitching Moment Coefficient")
    println("Cm: $Cm")
    # println("Yawing Moment Coefficient")
    println("Cn: $Cn")
    # println("Rolling Rate Coefficient")
    println("p̄: $p̄")
    # println("Pitching Rate Coefficient")
    println("q̄: $q̄")
    # println("Yawing Rate Coefficient")
    println("r̄: $r̄")
end

end