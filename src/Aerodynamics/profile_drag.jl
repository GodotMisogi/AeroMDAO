## Wetted-area method
#======================================================#

# Raymer form factor for wing, pylon, nacelle, fuselage? (Eq. 12.30)
form_factor_wing(x_c, t_c, Λ_m, M) = (1 + 0.6t_c / x_c + 100t_c^4) * (1.34M^0.18 * cos(Λ_m)^0.28)

function form_factor(wing :: HalfWing, M)
    # (x/c)_max, (t/c)_max, (t/c)_max sweep angles
    xcs, tcs, sweeps = max_tbyc_sweeps(wing, 60)
    Kf 				 = form_factor_wing.(xcs, tcs, sweeps, M)
end

# Schlichting averaged skin-friction coefficients
cf_lam(Re_c, k_lam = 1.) = 1.328 / √(Re_c * k_lam)
cf_turb(Re_c, M) = 0.455 / log10(Re_c)^2.58 / (1 + 0.144M^2)^0.65
cf_schlichting(Re, Re_xtr, M) = max(cf_lam(Re), cf_turb(Re, M) - (Re_xtr / 320 - 39) / Re)

function wetted_area_drag(wing :: HalfWing, x_tr, V, ρ, a_ref = 330., μ = 1.5e-5)
    # Chord processing
    chords      = wing.chords
    mean_chords = fwdsum(chords) / 2
    
    # Skin-friction coefficients
    M 			= V / a_ref
    Re_c 		= reynolds_number.(ρ, V, mean_chords, μ)
    Re_xtr 		= reynolds_number.(ρ, V, mean_chords .* x_tr, μ)
    cfs 		= cf_schlichting.(Re_c, Re_xtr, M)
    
    # Wetted areas
    S_wets 		= @. mean_chords * wing.spans / cos(wing.dihedrals)

    # Form factors
    K_fs 		= form_factor(wing, M)

    # Profile drag
    Dp_by_q	 	= sum(@. cfs * S_wets * K_fs)
end

profile_drag_coefficient(wing :: HalfWing, x_tr, V, rho_ref, a_ref, area_ref, μ) = wetted_area_drag(wing, x_tr, V, rho_ref, a_ref, μ) / area_ref
profile_drag_coefficient(wing :: Wing, x_tr, V, rho_ref, a_ref, area_ref, μ) = profile_drag_coefficient(wing.left, x_tr, V, rho_ref, a_ref, area_ref, μ) + profile_drag_coefficient(wing.right, x_tr, V, rho_ref, a_ref, area_ref, μ)

## Local-friction and local-dissipation method