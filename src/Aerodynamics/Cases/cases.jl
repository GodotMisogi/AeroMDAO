## 3D cases
#==========================================================================================#

"""
    print_info(wing :: Union{Wing, HalfWing})

Print the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
function print_info(wing :: Union{Wing, HalfWing}, head = ""; browser = false)
    labels = [ "Span (m)", "Area (m²)", "MAC (m)", "Aspect Ratio"]
    wing_info = [ info(wing)... ]
    data = [ labels wing_info ]
    header = [ head, "Value" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true)
    if browser
        pretty_table(String, data, header, alignment = [:c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, header, alignment = [:c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end


"""
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, freestream :: Freestream, r_ref, ρ = 1.225; symmetry = false)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, with an option for symmetry.
"""
function solve_case(horseshoe_panels :: Array{<: Panel3D}, normals, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref, r_ref)
    # Make horseshoes and collocation points
    horseshoes = horseshoe_line.(horseshoe_panels)
    collocation_points = horseshoe_point.(horseshoe_panels)

    # Solve system
    Γs = reshape(solve_system(horseshoes[:], collocation_points[:], normals[:], U, Ω), size(horseshoe_panels))

    # Compute forces and moments
    surface_forces, surface_moments, trefftz_force = case_dynamics(Γs, horseshoes, U, α, β, Ω, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(surface_forces, surface_moments, trefftz_force, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, Γs, horseshoes
end

"""
    solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref; span_num :: Integer = 5, chord_num :: Integer = 10)

Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
"""
function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, r_ref = [0.25, 0, 0], area_ref = 1, chord_ref = 1, span_ref = 1, mu_ref = 1.5e-5, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, viscous = false, a_ref = 330., x_tr = 0.3)
    # Unpack Freestream
    U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Determine spanwise panel distribution
    span_nums = number_of_spanwise_panels(wing, span_num)

    # Compute panels and normals
    horseshoe_panels, camber_panels = vlmesh_wing(wing, span_nums, chord_num)
    normals = panel_normal.(camber_panels)

    # Compute forces and moments
    nearfield_coeffs, farfield_coeffs, CFs, CMs, Γs, horseshoes = solve_case(horseshoe_panels, normals, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref, r_ref)
    
    # Viscous drag evaluation
    if viscous
        # Compute profile drag via wetted-area equivalent skin-friction method
        CDv = profile_drag_coefficient(wing, x_tr, freestream.V, rho_ref, a_ref, area_ref, mu_ref)

        # Add profile and total drag coefficients
        nf_coeffs = [ nearfield_coeffs[1] + CDv; CDv; nearfield_coeffs ]
        ff_coeffs = [  farfield_coeffs[1] + CDv; CDv; farfield_coeffs  ]
    else
        nf_coeffs = nearfield_coeffs
        ff_coeffs = farfield_coeffs
    end

    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs
end

function solve_case(components :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft", print = false, print_components = false) where T <: Real
    # Get panels
    meshes = values(components)
    
    # Flattening for VLM
    horseshoe_panels = first.(meshes)
    normals      	 = last.(meshes)
    horsies          = vcat(vec.(horseshoe_panels)...)
    normies          = vcat(vec.(normals)...)

    # Get required vortex lattice variables, i.e. horseshoes, collocation points and normals
    horseshoes = horseshoe_line.(horsies)
    collocation_points = horseshoe_point.(horsies)
    horseshoes_arr = [ horseshoe_line.(horses) for horses in horseshoe_panels ]
    
    # Unpack Freestream
    U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Solve system
    Γs = solve_system(horseshoes, collocation_points, normies, U, Ω)

    # Reshaping
    panel_sizes = size.(horseshoe_panels)
    panel_inds 	= [ 0; cumsum(prod.(panel_sizes)) ]
    Γs_arr 		= reshape_array(Γs, panel_inds, panel_sizes)

    # Compute forces and moments
    results = case_dynamics.(Γs_arr, horseshoes_arr, Ref(Γs), Ref(horseshoes), Ref(U), α, β, Ref(Ω), rho_ref, Ref(r_ref))

    # Components' non-dimensional forces and moments
    data = [ evaluate_coefficients(dyn..., U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref) for dyn in results ]

    nf_comp_coeffs = getindex.(data, 1)
    ff_comp_coeffs = getindex.(data, 2)
    CFs            = getindex.(data, 3)
    CMs            = getindex.(data, 4)

    # Aircraft's non-dimensional forces and moments
    nf_coeffs = reduce((x, y) -> x .+ y, nf_comp_coeffs) # Sum nearfield coefficients
    ff_coeffs = reduce((x, y) -> x .+ y, ff_comp_coeffs) # Sum farfield coefficients
    name_CFs  = vcat(vec.(CFs)...)						 # Collect surface force coefficients
    name_CMs  = vcat(vec.(CMs)...)						 # Collect surface moment coefficients

    # Printing
    if print_components; print_coefficients.(nf_comp_coeffs, ff_comp_coeffs, keys(components)) end
    if print; 			 print_coefficients(nf_coeffs, ff_coeffs, name) 					   end

    # Dictionary assembly
    name_data = (nf_coeffs, ff_coeffs, name_CFs, name_CMs, horsies, normies, horseshoes, Γs)
    comp_data = tuple.(nf_comp_coeffs, ff_comp_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes_arr, Γs_arr)

    names 	= [ name						 ; # Aircraft name
                (collect ∘ keys)(components) ] # Component names
    data  	= [ name_data ;	# Aircraft data
                comp_data ]	# Component data

    Dict(names .=> data)
end

streamlines(freestream :: Freestream, points, horseshoes, Γs, length, num_steps :: Integer) = VortexLattice.streamlines.(points, Ref(velocity(freestream)), Ref(freestream.omega), Ref(horseshoes), Ref(Γs), Ref(length), Ref(num_steps))