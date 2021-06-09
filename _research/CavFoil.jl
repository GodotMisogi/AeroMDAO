module CavFoil

include("DoubletSource.jl")

using .DoubletSource

function cavity_split_foil(foil :: Foil, cavity_start, cavity_end)
    upper, lower = split_foil(foil.coords)

    xs = [ first(coords) for coords in upper ]
    c = (abs ∘ map)(-, extrema(xs))
    cavity_le = cavity_start
    l = (cavity_end - cavity_start) / 100 / c

    leading_cav_upper = filter(x -> first(x) <= cavity_start * c / 100, upper)
    wetted_cav_upper = filter(x -> first(x) >= cavity_start * c / 100 && first(x) <= cavity_end * c / 100, upper)
    cav_trailing_upper = filter(x -> first(x) >= cavity_end * c / 100, upper)
    leading_cav_lower = filter(x -> first(x) <= cavity_start * c / 100, lower)
    wetted_cav_lower = filter(x -> first(x) >= cavity_start * c / 100 && first(x) <= cavity_end * c / 100, lower)
    cav_trailing_lower = filter(x -> first(x) >= cavity_end * c / 100, lower)

    return (leading_cav_upper, wetted_cav_upper, cav_trailing_upper, leading_cav_lower, wetted_cav_lower, cav_trailing_lower)
end

function cavity_length_foil(foil :: Foil, cavity_start, cavity_end, N_panels = 40, N_cavpanels = 10)

    leading_cav_upper, wetted_cav_upper, cav_trailing_upper, leading_cav_lower, wetted_cav_lower, cav_trailing_lower = cavity_split_foil(foil, cavity_start, cavity_end)

    wet_foil = cosine_foil(wetted_cav_lower, N_cavpanels)

    # Upper: Leading (cavity) to trailing
    wu = cosine_interp(wetted_cav_upper, Int(trunc(N_cavpanels)))
    tu = cosine_interp(cav_trailing_upper, Int(trunc((N_panels - N_cavpanels)/3)))
    # Lower: Leading (cavity) to trailing
    wl = cosine_interp(wetted_cav_lower, Int(trunc((N_panels - N_cavpanels)/3)))
    tl = cosine_interp(cav_trailing_lower, Int(trunc((N_panels - N_cavpanels)/3)))
    # Concatenate coordinates and reorder
    wet_upper = [wu; tu]
    wet_lower = [wl; tl]
end

# Cavity termination models 
cavity_pressure(sf, sl, A = 0.5, ν = 1.0, λ = 0.1) = sf < (1 - λ) * sl ? 0 : A * ((sf - (1 - λ) * sl) / (sl - (1 - λ) * sl))^ν

function cavity_influence_matrix(wetted_panels :: Vector{Panel2D}, wetted_cavpanels :: Vector{Panel2D}, cavity_panels :: Vector{Panel2D})
    # Influence on wetted panels
    doublets_wetwet   = doublet_matrix(wetted_panels, wetted_panels)
    sources_wetcav    = source_matrix(wetted_panels, cavity_panels)

    # Influence on cavity panels
    doublets_cavwet   = doublet_matrix(cavity_panels, wetted_panels)
    sources_cavcav    = source_matrix(cavity_panels, cavity_panels)

    # Kutta condition
    kutta             = kutta_condition(wetted_panels)
    wake_wet          = wake_vector(wetted_panels)
    wake_cav          = wake_vector(cavity_panels)

    # Cavity transition condition
    doublets_wetcav   = doublet_matrix(wetted_panels, cavity_panels)
    doublets_cavcav   = doublet_matrix(cavity_panels, cavity_panels)
  
    sf                = panel_length.(wetted_cavpanels)
    sl                = sum(panel_length.(cavity_panels))
    trans             = cumsum(panel_length.(cavity_panels)) .* (1 .- cavity_pressure.(cumsum(sf), sl))

    trans_wetcav      = doublets_wetcav * trans
    trans_cavcav      = doublets_cavcav * trans

    # Cavity closure condition
    closure_cav       = cavity_length ./ (1 - cavity_pressure.(cumsum(sf), sl))

    # Cavity leading edge extrapolation
    leading_cavpanel  = findfirst(panel -> panel === first(cavity_panels), panels)
    prev_terms        = [1,2,3]
    leading           = [ zeros(length(panels) - leading_cavpanel - length(prev_terms)); prev_terms; zeros(leading_cavpanel - 1) ]
    leading_dub       = sum(doublets_wetcav, dims = 2)
    leading_cav       = sum(doublets_cavcav, dims = 2)

    # Aerodynamic Influence Coefficient matrix
    AIC               = [ doublets_wetwet  sources_wetcav  trans_wetcav  leading_dub  wake_wet ;
                          doublets_cavwet  sources_cavcav  trans_cavcav  leading_cav  wake_cav ;
                          zeros(num_wet)'    closure_cav'        0            0           0    ;
                             leading'      zeros(num_cav)'       0            0           0    ;
                              kutta'       zeros(num_cav)'       0            0           0    ]

end

function cavity_boundary_vector(wetted_panels :: Vector{Panel2D}, wetted_cavpanels :: Vector{Panel2D}, cavity_panels :: Vector{Panel2D}, freestream :: Uniform2D)
    # Source terms
    boundary_sources_wet    = boundary_condition(wetted_panels, freestream)
    boundary_sources_cav    = - source_matrix(cavity_panels, wetted_panels) * source_strengths(wetted_panels, freestream)

    # Doublet terms
    bound_cav               = [ potential(freestream, colpoint...) - potential(freestream, (p1 ∘ first)(cavity_panels)...) for colpoint in collocation_point.(cavity_panels) ]
    
    boundary_doublets_wet  = doublet_matrix(cavity_panels, wetted_panels) * bound_cav
    boundary_doublets_cav  = doublet_matrix(cavity_panels, cavity_panels) * bound_cav

    # Cavity closure
    sf                      = panel_length.(wetted_cavpanels)
    sl                      = sum(panel_length.(cavity_panels))
    sources_cav             = source_strengths(cavity_panels, freestream)
    closure_cav             = cavity_length ./ (1 - cavity_pressure.(cumsum(sf), sl))
    source_closure          = - sum(sources_cav .* closure_cav)

    # Matrix assembly
    RHS                     = [ boundary_doublets_wet + boundary_sources_wet ;
                                boundary_doublets_cav + boundary_sources_cav ;
                                               source_closure                ;
                                                     0                       ;
                                                     0                       ]
end    


function solve_strengths(panels :: Vector{Panel2D}, cavity_panels :: Vector{Panel2D}, freestream :: Uniform2D)
    wetted_panels, wetted_cavpanels = cavity_split_panels(panels)

    # Matrix solution
    AIC       = cavity_influence_matrix(wetted_panels, wetted_cavpanels, cavity_panels)
    RHS       = cavity_boundary_vector(wetted_panels, wetted_cavpanels, cavity_panels, freestream)
    strengths = AIC \ RHS

    φs        = strengths[1:num_wetted]
    ∂φ∂ns     = strengths[num_wetted+1:end-2]

    # Cavity velocity, cavitation number and wake panel strength
    qc        = strengths[end - 1]
    σ         = qc^2 / freestream.V^2 - 1

    println("Cavity Velocity: ", qc, ", Cavitation Number: ", σ)
end

function cavity_height(cavity_panels :: Vector{Panel2D})
    heights = zeros(num_cavity)
    cavity_length = sum([ panel.length for panel in cavity_panels ])
    for n in 1:num_cavity
        sfj += wetted_cavpanels[n].length
        heights[n] = cavity_length * (sum(u .* cavity_panels[n].normal) + cavity_panels[n].source_strength) / (qc * (1 - cavPressure(sfj, sl)))
    end

    for (h1, h2, p1, p2) in zip(heights[1:end - 1], heights[2:end], cavity_panels[1:end - 1], cavity_panels[2:end])
        # Update coordinates
        angle = p1.normal .+ p2.normal
        p1.xe = p2.xs += (h1 + h2) / 2 * angle[1]
        p1.ye = p2.ys += (h1 + h2) / 2 * angle[2]

        # Panel 1 update
        p1.length = mag([p1.xe - p1.xs, p1.ye - p1.ys])
        p1.angle = atan(p1.ye - p1.ys, p1.xe - p1.xs)
        p1.normal = (-sin(p1.angle), cos(p1.angle))
        p1.tangent = (cos(p1.angle), sin(p1.angle))
        p1.xc, p1.yc = (p1.xe + p1.xs) / 2.0, (p1.ye + p1.ys) / 2.0

        # Panel 2 update
        p2.length = mag([p2.xe - p2.xs, p2.ye - p2.ys])
        p2.angle = atan(p2.ye - p2.ys, p2.xe - p2.xs)
        p2.normal = (-sin(p2.angle), cos(p2.angle))
        p2.tangent = (cos(p2.angle), sin(p2.angle))
        p2.xc, p2.yc = (p2.xe + p2.xs) / 2.0, (p2.ye + p2.ys) / 2.0
    end

    return heights
end

end