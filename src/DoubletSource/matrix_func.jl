"""
doublet_matrix(panels_1, panels_2)

Create the matrix of doublet potential influence coefficients between pairs of panels_1 and panels_2.
"""
doublet_matrix(panels_1 :: Vector{<: Panel2D}, panels_2 :: Vector{<: Panel2D}) = [ ifelse(panel_i === panel_j, 0.5, doublet_influence(panel_j, panel_i)) for panel_i in panels_1, panel_j in panels_2 ]

"""
	kutta_condition(panels)

Create the vector describing Morino's Kutta condition given Panel2Ds.
"""
kutta_condition(panels :: Vector{<: Panel2D}) = [ 1 -1 zeros(length(panels) - 4)' 1 -1 ]

"""
	wake_vector(wake_panel, panels)

Create the vector of doublet potential influence coefficients from the wake on the panels given the wake panel and the array of Panel2Ds.
"""
wake_vector(woke_panel :: Panel2D, panels :: Vector{<: Panel2D}) = doublet_influence.(Ref(woke_panel), panels)

"""
	influence_matrix(panels, wake_panel)

Assemble the Aerodynamic Influence Coefficient matrix consisting of the doublet matrix, wake vector, Kutta condition given Panel2Ds and the wake panel.
"""
influence_matrix(panels :: Vector{<: Panel2D}, woke_panel :: Panel2D) =
	[ doublet_matrix(panels, panels)  wake_vector(woke_panel, panels) ;
		kutta_condition(panels)						0				  ]


"""
	source_matrix(panels_1, panels_2)

Create the matrix of source potential influence coefficients between pairs of `panels_1` and `panels_2`.
"""
source_matrix(panels_1 :: Vector{<: Panel2D}, panels_2 :: Vector{<: Panel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

"""
	source_strengths(panels, freestream)

Create the vector of source strengths for the Dirichlet boundary condition ``\\sigma = \\vec U_{\\infty} \\cdot \\hat{n}`` given Panel2Ds and a Uniform2D.
"""
source_strengths(panels :: Vector{<: Panel2D}, u) = dot.(Ref(u), panel_normal.(panels))

"""
	boundary_vector(panels, u)

Create the vector for the boundary condition of the problem given an array of Panel2Ds and velocity ``u``.
"""
boundary_vector(panels :: Vector{<: Panel2D}, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

boundary_vector(colpoints, u) = [ dot.(colpoints, Ref(u)); 0 ]

function boundary_vector(panels :: Vector{<: Panel2D}, wakes :: Vector{<: Panel2D}, u) 
	source_panels = [ panels; wakes ]
	[ - source_matrix(panels, source_panels) * source_strengths(source_panels, u); 0 ]
end

"""
	solve_strengths(panels, u, wakes, bound)

Solve the system of equations ``[AIC][\\phi] = [\\vec{U} \\cdot \\hat{n}] - B[\\sigma]`` condition given the array of Panel2Ds, a velocity ``\\vec U``, a vector of wake `Panel2D`s, and an optional named bound for the length of the wake.
"""
function solve_strengths(panels :: Vector{<: Panel2D}, u, wakes :: Vector{<: Panel2D}; bound = 1e2) 
	AIC  = influence_matrix(panels, wake_panel(panels, bound))
	boco = boundary_vector(panels, wakes, u)

	AIC \ boco 
end

"""
	solve_strengths(panels, u, sources, bound)

Solve the system of equations ``[AIC][\\phi] = [\\vec{U} \\cdot \\hat{n}] - B[\\sigma]`` condition given the array of Panel2Ds, a velocity ``\\vec U``, a condition whether to disable source terms (``\\sigma = 0``), and an optional named bound for the length of the wake.
"""
function solve_strengths(panels :: Vector{<: Panel2D}, u, sources :: Bool; bound = 1e2) 
	AIC  = influence_matrix(panels, wake_panel(panels, bound))
	boco = boundary_vector(ifelse(sources, panels, collocation_point.(panels)), u)

	AIC \ boco 
end

"""
	lift_coefficient(panels, φs, u, sources :: Bool)

Compute the pressure coefficient and lift coefficient distributions given the array of `Panel2D`s, the associated doublet strengths ``\\phi``s, the velocity ``u``, and a condition whether to disable source terms (``\\sigma = 0``).
"""
function aerodynamic_coefficients(panels :: Vector{<: Panel2D}, φs, u, sources :: Bool)
	speed = norm(u)
	
	Δrs   = midpair_map(panel_dist, panels)
	Δφs   = -midpair_map(-, φs[1:end-1])
	vels  = ifelse(sources, panel_velocity.(Δφs, Δrs, Ref(u), panel_tangent.(panels)), Δφs ./ Δrs)
	cps   = @. pressure_coefficient(speed, vels)
	cls   = @. lift_coefficient(cps, Δrs / 2, panel_angle(panels))
	
	cps, cls
end