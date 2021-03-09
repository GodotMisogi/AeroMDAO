module NonDimensional

using LinearAlgebra: norm

export dynamic_pressure, force_coefficient, moment_coefficient, moment_coefficients, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics, reynolds_number

"""
	dynamic_pressure(ρ, V)

Compute the dynamic pressure given density ``\\rho`` and speed ``V``.
"""
dynamic_pressure(ρ, V) = 1/2 * ρ * V^2

"""
	force_coefficient(force, q, S)

Compute the non-dimensional force coefficient given force, dynamic pressure ``q``, and area ``S``.
"""
force_coefficient(force, q, S) = force / (q * S)

"""
	moment_coefficient(moment, q, S, ref_length)

Compute the non-dimensional moment coefficient given force, dynamic pressure ``q``, area ``S``, and reference length.
"""
moment_coefficient(moment, q, S, ref_length) = force_coefficient(moment, q, S) / ref_length

"""
	rate_coefficient(Ω, V, ref_length)

Compute the non-dimensional angular velocity coefficient given angular speed ``\\Omega``, reference speed ``V``, and reference length.
"""
rate_coefficient(Ω, V, ref_length) = Ω * ref_length / 2V

"""
	pressure_coefficient(force, ρ, V, S)

Compute the pressure coefficient given force, reference density ``\\rho``, speed ``V`` and area ``S``.
"""
pressure_coefficient(force, ρ, V, S) = force_coefficient(force, dynamic_pressure(ρ, V), S)

pressure_coefficient(mag, vels) = 1 - (norm(vels) / mag)^2

function moment_coefficient(moment, q, S, b, c)
	Cl = moment_coefficient(moment[1], q, S, b)
	Cm = moment_coefficient(moment[2], q, S, c)
	Cn = moment_coefficient(moment[3], q, S, b)

	Cl, Cm, Cn
end

function rate_coefficient(Ω, V, b, c)
	p̄	= rate_coefficient(Ω[1], V, b) 
	q̄	= rate_coefficient(Ω[2], V, c)
	r̄	= rate_coefficient(Ω[3], V, b)

	p̄, q̄, r̄
end
"""
	aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)

Prints the relevant aerodynamic/flight dynamics information given a net force, moment and angular rates ``\\Omega`` with reference speed ``V``, area ``S``, span ``b``, chord ``c``, and density ``\\rho``.
"""
function aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)
	q 			= dynamic_pressure(ρ, V)

	CDi, CY, CL = force_coefficient(force, q, S) 
	Cl, Cm, Cn	= moment_coefficient(moment, q, S, b, c)
	p̄, q̄, r̄	   = rate_coefficient(Ω, V, b, c)

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

reynolds_number(ρ, V, c, μ) = ρ * V * c / μ
reynolds_number(V, c, ν) = V * c / ν

end