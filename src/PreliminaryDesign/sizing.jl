# Generic
dynamic_pressure(ρ, V) = 1/2 * ρ * V^2
stall_speed(wing_loading, CL_max, ρ = 1.225) = 2 * wing_loading/(ρ * CL_max)

# Drag polar
span_efficiency_factor(e, AR) = 1 / (π * e * AR)
drag_polar(CD0, k, CL) = CD0 + k * CL^2

# Wing loading
wing_loading_stall_speed(V_stall, CL_max, ρ = 1.225) = 1/2 * ρ * V_stall^2 * CL_max


## Weight estimation and performance
#==========================================================================================#

# Power loading
power_loading(tbW, V, η_p) = tbW * V / η_p

flat_plate_cd(α) = 2(sin(α))^3

figure_of_merit(thrust) = 0.4742 * thrust^0.0793

hover_speed(T, S_p, ρ) = sqrt(1/2 * T / (ρ * S_p))

induced_velocity(climb_rate_vtol, v_h) = let x = climb_rate_vtol / (2v_h); -x + sqrt(x^2 + 1) end

disk_loading(M_TO) = 3.2261M_TO + 74.991

prop_area(W_TO, DL, η_prop) = W_TO / (DL * η_prop)