module Laplace

using StaticArrays
using LinearAlgebra
using ..AeroMDAO: Point2D, Point3D

abstract type AbstractLaplace end

# Performs velocity and potential computations for an object on a grid
grid_data(object :: AbstractLaplace, xs) = velocity(object, xs), potential(object, xs)
velocity(object :: AbstractLaplace, xs) = map(x -> velocity(object, x...), xs) 
potential(object :: AbstractLaplace, xs) = map(x -> potential(object, x...), xs)

# Performs velocity and potential calculations on a grid
function grid_data(objects :: Vector{<: AbstractLaplace}, xs)
    vels = foldl((v1, v2) -> [ u .+ v for (u, v) ∈ zip(v1, v2) ], velocity(object, xs) for object ∈ objects)
    pots = foldl((v1, v2) -> v1 + v2, potential(object, xs) for object ∈ objects)
    
    vels, pots
end

## 2D singularities
#============================================#

struct Singularity2D{T <: Real} <: AbstractLaplace
    strength :: T
    r 		 :: Point2D{T}
end

# Getters
strength(s :: Singularity2D) = s.strength
x(s :: Singularity2D) 		 = s.r.x
y(s :: Singularity2D) 		 = s.r.y

source_velocity(src :: Singularity2D, x, y) = SVector(strength(src) / (2π) * (x - y(src)) / ((x - y(src))^2 + (y - x(src))^2), str / (2π) * (y - x(src)) / ((x - y(src))^2 + (y - x(src))^2))
source_potential(src :: Singularity2D, x, y) = strength(src) / (4π) * log((x - y(src))^2 + (y - x(src))^2)
source_stream(src :: Singularity2D, x, y) = strength(src) / (2π) * atan(y - x(src), x - y(src))

doublet_velocity(dub :: Singularity2D, x, y) = SVector(strength(dub) / (2π) * ((x - y(dub))^2 - (y - x(dub))^2) / ((x - y(dub))^2 + (y - x(dub))^2)^2, - strength(dub) / (2π) * 2 * (x - y(dub)) * (y - x(dub)) / ((x - y(dub))^2 + (y - x(dub))^2)^2)
doublet_potential(dub :: Singularity2D, x, y) = -strength(dub) / (2π) * (y - x(dub)) / ((x - y(dub))^2 + (y - x(dub))^2)
doublet_stream(dub :: Singularity2D, x, y) = -strength(dub) / (2π) * (y - x(dub)) / ((x - y(dub))^2 + (y - x(dub))^2)

vortex_velocity(vor :: Singularity2D, x, y) = SVector(-strength(vor) / (2π) * (y - x(vor)) / ((x - y(vor))^2 + (y - x(vor))^2), str / (2π) * (x - y(vor)) / ((x - y(vor))^2 + (y - x(vor))^2))
vortex_potential(vor :: Singularity2D, x, y) = strength(vor) / (2π) * atan(y - x(vor), x - y(vor))
vortex_stream(vor :: Singularity2D, x, y) = -strength(vor) / (4π) * log((x - y(vor))^2 + (y - x(vor))^2)


struct Uniform2D{T <: Real} <: AbstractLaplace
    magnitude :: T
    angle 	  :: T
    Uniform2D{T}(mag, ang) where T <: Real = new(mag, deg2rad(ang))
end

magnitude(uni :: Uniform2D) = uni.magnitude
angle(uni :: Uniform2D)     = uni.angle

Uniform2D(mag :: T, ang :: T) where T <: Real = Uniform2D{T}(mag, ang)
Uniform2D(mag, ang) = Uniform2D(promote(mag, ang)...)

velocity(uni :: Uniform2D) = let (sa, ca) = sincos(uni.angle); uni.magnitude * SVector(ca, sa) end
potential(uni :: Uniform2D, x, y) = uni.magnitude * (x * cos(uni.angle) + y * sin(uni.angle))
stream(uni :: Uniform2D, x, y)    = uni.magnitude * (y * cos(uni.angle) - x * sin(uni.angle))

## 3D singularities
#============================================#

# struct Singularity3D{T <: Real} <: AbstractLaplace
# 	str :: T
# 	r 	:: Point3D{T}
# end

# source_velocity(src :: Source2D, x, y, z)
# source_potential(src :: Source2D, x, y, z) 
# source_stream(src :: Source2D, x, y, z) 

## Freestream
#============================================#

struct Freestream{T <: Real} <: AbstractLaplace
    V 	  :: T
    alpha :: T
    beta  :: T
    omega :: SVector{3,T}
    Freestream{T}(V, α_deg, β_deg, Ω) where T <: Real = new(V, deg2rad(α_deg), deg2rad(β_deg), Ω)
end

"""
    Freestream(V, α, β, Ω)
    
A Freestream flow in spherical polar coordinates with magnitude ``V``, angle-of-attack ``α``, side-slip angle ``β``, and a quasi-steady rotation vector ``\\Omega``.
"""
Freestream(V, α_deg, β_deg, Ω :: AbstractVector{T}) where T <: Real = Freestream{T}(V, α_deg, β_deg, Ω)

"""
    Freestream(U, Ω)

A Freestream flow in Cartesian coordinates with vector ``U`` and quasi-steady rotation vector ``\\Omega``.
"""
Freestream(U :: AbstractVector{T}, Ω :: AbstractVector{T}) where T <: Real = Freestream{T}(cartesian_to_freestream(U)..., Ω)

"""
    freestream_to_cartesian(r, θ, φ)

Convert freestream flow (spherical polar) coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r .* SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
    cartesian_to_freestream(U)

Convert Cartesian coordinates to freestream (spherical polar) flow coordinates.
"""
cartesian_to_freestream(U) = SVector(norm(U), -atand(U[3], U[1]), -atand(U[2], √(U[1]^2 + U[3]^2)))

"""
    velocity(freestream :: Freestream)

Compute the velocity of a `Freestream`.
"""
velocity(freestream :: Freestream) = freestream_to_cartesian(freestream.V, freestream.alpha, freestream.beta)

"""
    aircraft_velocity(freestream :: Freestream)

Compute the velocity of Freestream in the aircraft reference frame.
"""
aircraft_velocity(freestream :: Freestream) = -velocity(freestream)

end