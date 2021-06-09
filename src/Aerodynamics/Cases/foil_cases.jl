## 2D cases
#==========================================================================================#

"""
    solve_case(foil :: Foil, uniform :: Uniform2D; sources = false, wake_length = 1e3, num_panels :: Integer = 60)

Evaluate a doublet-source case given a `Foil` with a `Uniform2D`, with optional named arguments to specify whether the source terms are non-zero, the length of the wake, and the number of panels for the analysis.
"""
function solve_case(foil :: Foil, freestream :: Uniform2D; viscous = false, sources = false, wake_length = 1e5, num_panels :: Integer = 60, wake_panels = 15)
    panels = paneller(foil, num_panels)
    cls, cms, cps, cl_wake = solve_problem(panels, velocity(freestream), freestream.angle, ifelse(viscous, wake_panels, sources), wake_length)
    # cdp, cl = sincos(freestream.angle) .* cl_wake
    # cm      = cl * r_ref
    # coeffs    = (cdp, cl, cm)

    cl_wake, cls, cms, cps, panels
end