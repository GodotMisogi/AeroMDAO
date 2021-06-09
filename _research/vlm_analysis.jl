## VLM case setup
vlm = VLMAnalysis(aircraft, fs; 
                  rho_ref   = ρ, 
                  r_ref     = ref, 
                  area_ref  = S, 
                  span_ref  = b, 
                  chord_ref = c, 
                  name      = ac_name);


function vlm_analysis(aircraft, fs, ρ, ref, S, b, c)
    # Evaluate case
    data = 
    solve_case(aircraft, fs; 
               rho_ref   = ρ, 
               r_ref     = ref,
               area_ref  = S,
               span_ref  = b,
               chord_ref = c,
            #    print 	 = true
              );
        
    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = data["Aircraft"]
               
    SVector{3, Float64}(ff_coeffs[1:3]), SVector{3, Float64}(nf_coeffs[4:6])
end