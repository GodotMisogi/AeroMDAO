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
	v2 = @. v1 + p2 - p1
	v3 = @. v4 + p3 - p4
	
	v1, v2, v3, v4
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
Sums the velocities evaluated at a point `r` of vortex lines with constant strength Γ.
"""
sum_vortices(r, vortex_lines :: Array{<: Line}, Γ) = sum(velocity(r, line, Γ) for line ∈ vortex_lines)

"""
Computes the induced velocities at a point `r` of a Vortex Ring with constant strength Γ.
"""
velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)