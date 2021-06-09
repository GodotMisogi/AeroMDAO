module NonDimensional

using LinearAlgebra: norm
using PrettyTables

export dynamic_pressure, force_coefficient, moment_coefficient, moment_coefficients, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_coefficients, reynolds_number, print_derivatives

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
    moment_coefficient(moment, q, S, L)

Compute the non-dimensional moment coefficient given force, dynamic pressure ``q``, area ``S``, and reference length ``L``.
"""
moment_coefficient(moment, q, S, L) = force_coefficient(moment, q, S) / L

"""
    rate_coefficient(Ω, V, L)

Compute the non-dimensional angular velocity coefficient given angular speed ``\\Omega``, reference speed ``V``, and reference length ``L``.
"""
rate_coefficient(Ω, V, L) = Ω * L / 2V

"""
    pressure_coefficient(force, ρ, V, S)

Compute the pressure coefficient given force, reference density ``\\rho``, speed ``V`` and area ``S``.
"""
pressure_coefficient(force, ρ, V, S) = force_coefficient(force, dynamic_pressure(ρ, V), S)

pressure_coefficient(mag, vels) = 1 - (norm(vels) / mag)^2

moment_coefficient(moment, q, S, b, c) = moment_coefficient.(moment, q, S, [b, c, b])

rate_coefficient(Ω, V, b, c) = rate_coefficient.(Ω, V, [b, c, b])


"""
    aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)

Compute the relevant aerodynamic coefficients given a net force, moment and angular rates ``\\Omega`` with reference speed ``V``, area ``S``, span ``b``, chord ``c``, and density ``\\rho``.
"""
function aerodynamic_coefficients(force, moment, Ω, V, S, b, c, ρ)
    q = dynamic_pressure(ρ, V)

    [ force_coefficient(force, q, S) 		 ;
      moment_coefficient(moment, q, S, b, c) ;
      rate_coefficient(Ω, V, b, c) 			 ]
end


function print_coefficients(CL, CDi, CY, Cl, Cm, Cn, p̄, q̄, r̄)
    # println("Lift Coefficient")
    println("CL: $CL")
    # println("Drag Coefficient")
    println("CDi: $CDi")
    # println("Side Force Coefficient")
    println("CY: $CY")
    # println("Lift-to-Drag Ratio")
    println("L/D: $(CL ./ CDi)")
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

function print_coefficients(nf_coeffs, ff_coeffs, name = ""; browser = false)
    coeffs = [ ifelse(length(nf_coeffs) == 11, ["CD", "CDp"], []); [ "CDi", "CY", "CL", "Cl", "Cm", "Cn", "p̄", "q̄", "r̄" ] ]
    data = [ coeffs nf_coeffs [ ff_coeffs; fill("—", 6) ] ]
    head = [ name, "Nearfield", "Farfield" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true)
    if browser
        pretty_table(String, data, head, alignment = [:c, :c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, head, alignment = [:c, :c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end

function print_derivatives(derivs, name = ""; browser = false)
    coeffs = ["∂CD", "∂CY", "∂CL", "∂Cl", "∂Cm", "∂Cn"]
    nf_vars = [ "$name" "" "Nearfield" "Stability" "Derivatives" "" ; "" "∂α, 1/rad" "∂β, 1/rad" "∂p̄" "∂q̄" "∂r̄" ]
    nf_rows = [ coeffs derivs ]

    if browser
        pretty_table(String, nf_rows, nf_vars, alignment = :c, tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(color = "blue", font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(nf_rows, nf_vars, alignment = :c, tf = tf_compact, header_crayon = Crayon(bold = true), subheader_crayon = Crayon(foreground = :yellow, bold = true), highlighters = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true), vlines = :none, formatters = ft_round(8))
    end
end

reynolds_number(ρ, V, c, μ) = ρ * V * c / μ
reynolds_number(V, c, ν) = V * c / ν

end