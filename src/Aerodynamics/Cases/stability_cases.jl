## Stability derivative cases
function scale_inputs(freestream :: Freestream, b, c)
    # Building input vector
    x0 = [ freestream.alpha; 
           freestream.beta; 
           rate_coefficient(freestream.omega, freestream.V, b, c) ]

    # Unscaling non-dimensional rate coefficients
    scale = 2freestream.V .* [1/b, 1/c, 1/b]

    x0, scale
end

function print_case(nf, ff, derivs, comp)
    print_coefficients(nf, ff, comp)
    print_derivatives(derivs, comp)
end

function print_case(data, comp)
    nf, ff, derivs = data[comp]
    print_case(nf, ff, derivs, comp)
end

function solve_stability_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, name = "Wing", viscous = false, x_tr = 0.3, print = false)
    # Reference values and scaling inputs
    S, b, c = area_ref, span_ref, chord_ref
    x, scale = scale_inputs(freestream, b, c)

    # Closure to generate results with input vector
    function stab(x)
        fs 	 = Freestream(freestream.V, rad2deg(x[1]), rad2deg(x[2]), x[3:end] .* scale)
        data = solve_case(wing, fs;
                          rho_ref   = rho_ref, 
                          r_ref     = r_ref, 
                          area_ref  = S, 
                          span_ref  = b, 
                          chord_ref = c, 
                          chord_num = chord_num,
                          span_num  = span_num,
                          viscous 	= viscous,
                          x_tr 		= x_tr)

        [ data[1]; data[2] ] # Concatenate nearfield and farfield coefficients for DiffResults value
    end

    # Set up Jacobian system and evaluate
    y = ifelse(viscous, zeros(16), zeros(12))
    result = DiffResults.JacobianResult(y, x)
    result = ForwardDiff.jacobian!(result, stab, x)
    vars   = DiffResults.value(result)
    derivs = DiffResults.jacobian(result)

    # Gather results
    nf 	= ifelse(viscous, vars[1:11], vars[1:9]) 	  	# Nearfield coefficients
    ff 	= ifelse(viscous, vars[12:end], vars[10:end]) 	# Farfield coefficients
    dvs = ifelse(viscous, derivs[[1,4:8...],:], derivs[1:6,:])	# Nearfield stability derivatives, uses total drag for viscous cases just in case of generalisations

    # Print if needed
    if print; print_case(nf, ff, derivs, name) end

    nf, ff, dvs
end

function solve_stability_case(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft", print = false, print_components = false) where T <: Real
    # Reference values and scaling inputs
    S, b, c = area_ref, span_ref, chord_ref
    x, scale = scale_inputs(freestream, b , c)

    # Closure to generate results with input vector
    function stab(x)
        fs 	 = Freestream(freestream.V, rad2deg(x[1]), rad2deg(x[2]), x[3:end] .* scale)
        data = solve_case(aircraft, fs,
                          rho_ref   = rho_ref, 
                          r_ref     = r_ref, 
                          area_ref  = S, 
                          span_ref  = b, 
                          chord_ref = c,
                          name 		= name)

        # Creates array of nearfield and farfield coefficients for each component as a row vector.
        coeffs = hcat((let comps = data[name]; [ comps[1]; comps[2] ] end for name in names)...)
    end

    names 	= [ keys(aircraft)..., name ]
    num_comps = length(names)

    y 		= zeros(12, num_comps)
    result 	= DiffResults.JacobianResult(y, x)
    result 	= ForwardDiff.jacobian!(result, stab, x)

    vars 	= DiffResults.value(result)
    derivs 	= DiffResults.jacobian(result)

    # Dictionary assembly
    ranges  = (1:num_comps) .* 12                     # Painful hacking
    bounds  = zip([ 1; ranges[1:end-1] .+ 1], ranges) # Painful hacking
    data 	= Dict(name => (vars[1:9, i], vars[10:end, i], derivs[first(inds):last(inds)-6,:]) for (i, (name, inds)) in (enumerate ∘ zip)(names, bounds)) 
    
    # Printing
    if print;            print_case(data, name)       end
    if print_components; print_case.(Ref(data), names)end

    data
end


