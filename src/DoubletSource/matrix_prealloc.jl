"""
	kutta_condition!(AIC)

Set the Morino Kutta condition given the pre-allocated Aerodynamic Influence Coefficient Matrix `AIC`.
"""
function kutta_condition!(AIC)
	AIC[end,1] = 1
	AIC[end,2] = -1
	AIC[end,end-2] = 1
	AIC[end,end-1] = -1
	nothing
end

"""
	matrix_assembly!(AIC, boco, panels)

Build the AIC matrix and boundary condition vectors given the pre-allocated matrix `AIC`, pre-allocated boundary condition vector `boco`, and the array of `Panel2D`s.
"""
function matrix_assembly!(AIC, boco, panels, woke_panel :: Panel2D, u)
	for i in eachindex(panels)
		for j in eachindex(panels)
			boco[i] += boundary_condition(panels[j], panels[i], u)
			AIC[i,j] = ifelse(i == j, 0.5, doublet_influence(panels[j], panels[i]))
		end
		AIC[i,end] = doublet_influence(woke_panel, panels[i])
	end
	kutta_condition!(AIC)
	nothing
end

"""
	solve_strengths_prealloc(panels, u, bound = 1e3)

Solve the system of equations ``[AIC][\\phi] = [\\hat{U} \\cdot \\vec{n}]`` given the array of Panel2Ds, a freestream velocity vector ``u``, and an optional bound for the length of the wake. 
"""
function solve_strengths_prealloc(panels, u, bound = 1e3)
	n = length(panels) + 1
	AIC = zeros(n,n)
	boco = zeros(n)
	matrix_assembly!(AIC, boco, panels, wake_panel(panels, bound), u)

	AIC \ boco
end

"""
	panel_distances!(vec, panels, φs, u)

Compute the distances between adjacent panels given the pre-allocated `vec`, and an array of `Panel2D`s.
"""
function panel_distances!(vec, panels)
	vec[1] = panel_dist(panels[1], panels[2])
	vec[end] = panel_dist(panels[end-1], panels[end])
	for i in 2:length(panels) - 1
		vec[i] = panel_dist(panels[i-1], panels[i+1])
	end
	nothing
end

"""
	panel_velocities!(vec, panels, φs, u)

Compute the panel velocities given the pre-allocated `vec`, an array of `Panel2D`s, their associated doublet strengths ``\\phi``s, and a freestream velocity vector ``u``.
"""
function panel_velocities!(vec, panels, φs, u)
	vec[1] = panel_velocity(φs[2] - φs[1], panel_dist(panels[1], panels[2]), u, panel_tangent(panels[1]))
	vec[end] = panel_velocity(φs[end] - φs[end-1], panel_dist(panels[end], panels[end-1]), u, panel_tangent(panels[end]))
	for i in 2:length(panels) - 1
		vec[i] = panel_velocity(φs[i+1] - φs[i-1], panel_dist(panels[i+1], panels[i-1]), u, panel_tangent(panels[i]))
	end
	nothing
end

"""
	lift_coefficient_prealloc(panels, φs, u)

Compute the lift coefficient using pre-allocation, given `Panel2D`s, their associated doublet strengths ``\\phi``s, and a freestream velocity vector ``u``.
"""
function lift_coefficient_prealloc(panels, φs, u)
	panel_dists = (zeros ∘ length)(panels)
	panel_vels  = (zeros ∘ length)(panels)
	speed       = norm(u)
	
	panel_distances!(panel_dists, panels)
	panel_velocities!(panel_vels, panels, φs[1:end-1], u)
	cps = @. pressure_coefficient(speed, panel_vels)
	cls = @. lift_coefficient(cps, panel_dists / 2, panel_angle(panels))

	cl = sum(cls)
end