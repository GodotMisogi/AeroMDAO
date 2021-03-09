# Axis transformations
#==========================================================================================#

"""
	body_to_stability_axes(coords, α)

Converts coordinates into stability axes with angle ``\\alpha``.
"""
body_to_stability_axes(coords, α :: T) where T <: Real = RotY{T}(α) * coords
											
"""
	body_to_wind_axes(coords, α, β)

Converts coordinates from body axes to wind axes with angles ``\\alpha,~ \\beta``.
"""
body_to_wind_axes(coords, α :: T, β :: T) where T <: Real = RotZY{T}(β, α) * coords

"""
	body_to_wind_axes(coords, α, β)

Converts coordinates from wind axes to body axes with angles ``\\alpha,~ \\beta``.
"""
wind_to_body_axes(coords, α :: T, β :: T) where T <: Real = RotZY{T}(-α, -β) * coords

"""
	body_to_wind_axes(coords, α, β)

Converts coordinates from body axes to wind axes with angles ``\\alpha,~ \\beta``.
"""
body_to_wind_axes(vector, freestream :: Freestream) = body_to_wind_axes(vector, freestream.α, freestream.β)

"""
	body_to_wind_axes(coords, freestream :: Freestream)

Converts coordinates from wind axes to body axes given a vector and a `Freestream`. 
"""
wind_to_body_axes(vector, freestream :: Freestream) = wind_to_body_axes(vector, freestream.α, freestream.β)

"""
	reflect_xz(vector)

Reflects the y-coordinate of a given 3-dimensional vector about the x-z plane.
"""
reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])

"""
	project_yz(vector)

Projects the y-coordinate of a given 3-dimensional vector about the x-z plane.
"""
project_yz(vector) = SVector(0, vector[2], vector[3])

"""
	reflect_yz(line)

Reflects a Line onto the ``x``-``z`` plane of its coordinate system.
"""
reflect_xz(line :: Line) = Line((reflect_xz ∘ point2)(line), (reflect_xz ∘ point1)(line))

"""
	reflect_yz(horseshoe)

Reflects a Horseshoe onto the ``x``-``z`` plane of its coordinate system.
"""
reflect_xz(horseshoe :: Horseshoe) = (Horseshoe ∘ reflect_xz ∘ bound_leg)(horseshoe)

"""
	reflect_yz(panel :: Panel3D)

Reflects a Horseshoe onto the ``x``-``z`` plane of its coordinate system.
"""
reflect_xz(panel :: Panel3D) = Panel3D((reflect_xz ∘ point1)(panel), (reflect_xz ∘ point2)(panel), (reflect_xz ∘ point3)(panel), (reflect_xz ∘ point4)(panel))

"""
	project_yz(line)

Projects a Line onto the ``y``-``z`` plane of its coordinate system.
"""
project_yz(line :: Line) = Line((project_yz ∘ point1)(line), (project_yz ∘ point2)(line))

"""
	stability_flip(vector)

Reflects the x- and z- coordinates of a given 3-dimensional vector about the y-z and x-y planes respectively.
"""
stability_flip(vector) = SVector(-vector[1], vector[2], -vector[3])

"""
	body_to_wind_axes(line, freestream)

Transforms a Line from body to wind axes in a given `Freestream`.
"""
body_to_wind_axes(line :: Line, freestream :: Freestream) = Line(body_to_wind_axes(line.r1, freestream.α, freestream.β), body_to_wind_axes(line.r2, freestream.α, freestream.β)) 

"""
	body_to_stability_axes(force, moment, freestream)

Transforms forces and moments from body to stability axes in a given `Freestream`.
"""
body_to_stability_axes(force, moment, freestream :: Freestream) = body_to_stability_axes(force, freestream.α), body_to_stability_axes(stability_flip(moment), freestream.α), body_to_stability_axes(stability_flip(freestream.Ω), freestream.α)

"""
	body_to_wind_axes(force, moment, freestream)

Transforms forces and moments from body to wind axes in a given `Freestream`.
"""
body_to_wind_axes(force, moment, freestream :: Freestream) = body_to_wind_axes(force, freestream.α, freestream.β), body_to_wind_axes(stability_flip(moment), freestream.α, freestream.β), body_to_wind_axes(stability_flip(freestream.Ω), freestream.α, freestream.β)