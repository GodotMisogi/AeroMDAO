### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8771f120-8253-11eb-2e91-d70d9914670e
begin
	using Revise
	using LinearAlgebra
	using StaticArrays
	using Rotations
	using AeroMDAO
	using Plots
end

# ╔═╡ c39fe410-82fb-11eb-25ed-95925ebc31ba
using PlutoUI

# ╔═╡ 186b5640-82f8-11eb-1aee-e13872baed6f
md"# AeroMDAO - VLM Aircraft Analysis"

# ╔═╡ c55c41e0-827e-11eb-0380-d9d33c6e7148
md"## Aircraft Components"

# ╔═╡ d504db00-8311-11eb-2b50-5d6cc93e92a1
md"Airfoils"

# ╔═╡ c26fc67e-8311-11eb-07d3-05a06845d290
begin
	foil_path = "..\\data\\airfoil_database\\sd7037.dat"
	sd7037 = read_foil(foil_path)
end;

# ╔═╡ ad3d8a40-8253-11eb-3709-870c8ba3b76f
md"Wing"

# ╔═╡ ad4663e0-8253-11eb-0e02-31d856b9d844
begin
	wing_foils = Foil.([ sd7037, sd7037 ])
	wing_right = HalfWing(wing_foils,
						  [1.0, 0.6],
						  [0.0, 2.0],
						  [5.0],
						  [11.3],
						  [2.29]);
	wing 	 = Wing(wing_right, wing_right)
end;

# ╔═╡ f5bb8c76-9068-4065-9f73-edfa80778e95
HTML(print_info(wing, "Wing"; browser = true))

# ╔═╡ 8f953ff7-5d2e-49cf-9abb-697c397c89c6
begin
	wing_mac = mean_aerodynamic_center(wing)
	wing_pos = [0., 0., 0.]
end;

# ╔═╡ fc1aaaf8-4ec1-485f-ae16-fc21e87212a8
wing_plan = plot_wing(wing; position = wing_pos);

# ╔═╡ ad769b9e-8253-11eb-1946-f7f19dded706
md"Horizontal Tail"

# ╔═╡ ad80d4d0-8253-11eb-36cf-5373bb8741ed
begin
	htail_foil = Foil(naca4((0,0,1,2)))
	htail_right = HalfWing(fill(htail_foil, 2),
						   [0.7, 0.42],
						   [0.0, 0.0],
						   [1.25],
						   [0.],
						   [6.39])
	htail = Wing(htail_right, htail_right)

end;

# ╔═╡ 91048f97-1725-4266-8adf-18592130cdd2
HTML(print_info(htail, "Horizontal Tail"; browser = true))

# ╔═╡ befe81cf-8965-442e-b641-22627770d248
begin
	htail_mac	= mean_aerodynamic_center(htail)
	htail_pos	= [5., 0., 0.]
	α_h_i		= 0.
end;

# ╔═╡ 9bcf96fe-e5ef-4d26-bdbf-740ac883b82e
htail_plan = plot_wing(htail; position = htail_pos);

# ╔═╡ ada43b50-8253-11eb-2456-87a605b4c2da
md"Vertical Tail"

# ╔═╡ adc3f852-8253-11eb-3524-c1b3a12dc1e4
begin
	vtail_foil = Foil(naca4((0,0,0,9)))
	vtail = HalfWing(fill(vtail_foil, 2), 
					 [0.7, 0.42],
					 [0.0, 0.0],
					 [1.0],
					 [0.],
					 [6.39])
end;

# ╔═╡ c65f00f0-7596-483d-8fc0-6f0c4007248f
HTML(print_info(vtail, "Vertical Tail"; browser = true))

# ╔═╡ 7e83479d-e096-4315-a33e-3299f211b045
begin
	vtail_mac	= mean_aerodynamic_center(vtail) # NEEDS FIXING FOR ROTATION
	vtail_pos	= [5., 0., 0.]
end;

# ╔═╡ 509f54a9-ef72-40eb-8848-e1467abb5610
vtail_plan = plot_wing(vtail; position = vtail_pos);

# ╔═╡ 2d7f9d77-0c90-40b4-87de-74dfe036d7f8
md"## Static Stability"

# ╔═╡ 2ce19aa2-a158-4724-8c31-9e725096266c
begin
	S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
	x_w = wing_pos + [ wing_mac[1], 0, 0 ]
end

# ╔═╡ cac0286b-1535-436a-bab5-f779cae32998
md"Horizontal Tail Volume Coefficient"

# ╔═╡ 8d017053-31e6-4025-8047-d704e394db08
begin
	S_h = projected_area(htail)
	x_h = htail_pos + [ htail_mac[1], 0, 0 ]
	l_h = norm(x_h - x_w)
	V_h = S_h * l_h / (S * c)
end

# ╔═╡ 924ccf0c-eb30-4063-8856-d8f36b275ad9
md"Vertical Tail Volume Coefficient"

# ╔═╡ 1d37024b-3808-4911-9135-5fa53fa4133a
begin
	S_v = projected_area(vtail)
	x_v = vtail_pos + [ vtail_mac[1], 0, 0 ]
	l_v = norm(x_v - x_w)
	V_v = S_v * l_v / (S * b)
end

# ╔═╡ f1669dee-8254-11eb-3788-393adbf558ec
md"## Vortex Lattice Method Analysis"

# ╔═╡ d9e1dbc2-827e-11eb-160f-038db46faaa6
md"### Assembly"

# ╔═╡ add8b8d0-8253-11eb-0cdc-8f1e28af9e2c
begin
	wing_panels  = 	panel_wing(wing, [12], 6;
							   position = wing_pos
							  )
	htail_panels =	panel_wing(htail, [8], 4;
							   position	= htail_pos,
							   angle 	= deg2rad(α_h_i),
							   axis 	= [0., 1., 0.]
							  )
	vtail_panels = 	panel_wing(vtail, [8], 4;
							   position = vtail_pos,
							   angle    = π/2
							  )
end;

# ╔═╡ b79f2b23-1340-4e58-8adb-d3dd65fbc2c7
aircraft = Dict("Wing" 			  	=> wing_panels,
				"Horizontal Tail" 	=> htail_panels,
				"Vertical Tail"     => vtail_panels);

# ╔═╡ adf2a972-8253-11eb-05f7-0bc93fd90563
begin
	ρ 			= 1.225
	ref 		= x_w
	V, α, β 	= 10.0, 2.0, 0.0
	Ω 			= [0.0, 0.0, 0.0]
	freestream 	= Freestream(V, α, β, Ω)
end;

# ╔═╡ 66adab50-82f8-11eb-0603-f3a50a03c3b4
data = solve_case(aircraft, 
				   freestream; 
				   rho_ref     = ρ, 
				   r_ref       = ref, 
				   area_ref    = S, 
				   span_ref    = b, 
				   chord_ref   = c, 
				   name        = "My Aircraft"
				  );

# ╔═╡ 5b1ac8a0-8301-11eb-0e9f-7b4b8bb2b0bd
names = (collect ∘ keys)(data)

# ╔═╡ fc17a380-82f8-11eb-30b6-555448b1d939
comp = names[2]

# ╔═╡ 2b3361a0-82f8-11eb-2cf1-a5bfeb3b83f5
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = data[comp];

# ╔═╡ ab566bc0-82f8-11eb-2b30-b9d62de30027
HTML(print_coefficients(comp, nf_coeffs, ff_coeffs, browser = true))

# ╔═╡ c0d4b620-82fb-11eb-3735-935581c41b80
begin
	φ_s = @bind φ Slider(0:1e-2:90, default = 30)
	ψ_s = @bind ψ Slider(0:1e-2:90, default = 60)
	z_s = @bind z_limit Slider(0:1e-2:10, default = 1.)
	stream = @bind stream CheckBox()
	plotting = @bind plotting CheckBox()
	md"""
	Horizontal: $(φ_s)
	Vertical: $(ψ_s)
	
	*z*-scale: $(z_s)
	Plot: $plotting
	Streamlines: $stream
	"""
end

# ╔═╡ 35f40af0-8315-11eb-2821-efd5d319ffcb
gr(dpi = 300)

# ╔═╡ af888bb0-8253-11eb-2dd5-5d0739fd1f15
begin
	horseshoe_coords = plot_panels(horseshoe_panels[:])
	camber_coords    = plot_panels(camber_panels[:])
end;

# ╔═╡ e2fe51e0-8312-11eb-38e0-0d996d94ced1
md"## Trefftz Plane"

# ╔═╡ fdd49a2c-5c03-4d2a-9910-7aa3f4a88fc4
sum(CFs)

# ╔═╡ 19387ec0-834f-11eb-347f-69385fa70bd3
function trefftz_plane_forces(horseshoe_panels, CFs, ρ, V, α, β)
	# qs = 0.5 * ρ * V^2 * panel_area.(horseshoe_panels)
	trans_CFs = body_to_wind_axes.(CFs, deg2rad(α), deg2rad(β))
	sec_ys = getindex.(midpoint.(horseshoe_panels[1,:]), 2)
	CDis = sum(getindex.(trans_CFs, 1), dims = 1)[:]
	CYs	 = sum(getindex.(trans_CFs, 2), dims = 1)[:]
	CLs  = sum(getindex.(trans_CFs, 3), dims = 1)[:]
	
	sec_ys, CLs, CDis, CYs
end

# ╔═╡ e845a360-8312-11eb-1551-b3e2ecf7cbd2
begin
	sec_ys, CLs, CDis, CYs = trefftz_plane_forces(horseshoe_panels, CFs, ρ, V, α, β)
	plot1 = plot(sec_ys, CDis, marker = :dot, label = nothing, ylabel = "CDi")
	plot2 = plot(sec_ys, CYs, marker = :dot, label = nothing, ylabel = "CY")
	plot3 = plot(sec_ys, 2sum(Γs, dims = 1)[:] / (V * c), marker = :dot, label = nothing, xlabel = "y", ylabel = "CL Loading")
	plot(plot1, plot2, plot3, layout = (3, 1))
	plot!()
end

# ╔═╡ a49f0c90-8255-11eb-03da-e97868878c21
md"## Streamlines"

# ╔═╡ f3d5f660-82fb-11eb-156f-b9443de5879c
md"Seed 1"

# ╔═╡ fb8df9c0-82fb-11eb-0c21-1fb0470f0842
begin
	num_points = 50
	max_z = 0.1
	y = span(wing) / 2 - 0.05
	seed_1 = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))
end;

# ╔═╡ ffd1bdf0-82fb-11eb-142d-9bc50fafd5ed
md"Seed 2"

# ╔═╡ af01f690-8253-11eb-232e-bb504f4c3535
begin
	span_points = 50
	init        = leading_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points)
	dx, dy, dz  = 0, 0, 1e-3
	seed_2        = [ init .+ Ref([dx, dy, dz])  ; 
					  init .+ Ref([dx, dy, -dz]) ];
	distance = 5
	num_stream_points = 200
end;

# ╔═╡ f1e63450-82fb-11eb-305d-9bec0e0cd50a
streams = ifelse(stream, plot_streams(freestream, seed_2, horseshoes, Γs, distance, num_stream_points), nothing)

# ╔═╡ b04f5e20-8253-11eb-1263-934f17649375
begin
	if plotting 
		plot(xaxis = "x", yaxis = "y", zaxis = "z",
			 aspect_ratio = 1, 
			 camera = (φ, ψ),
			 zlim = (-0.1, z_limit),
			 )
		plot!.(camber_coords, color = :black, label = :none)
		if stream
			plot!.(streams, color = :green, label = :none)
		end
		plot!()
	end
end

# ╔═╡ bdd5da00-8277-11eb-08d7-0defcf2c7ba0
#= html"""<style>
main {
    max-width: 60%;
    align-self: flex-start;
    margin-left: 20%;
}
""" =#

# ╔═╡ Cell order:
# ╟─186b5640-82f8-11eb-1aee-e13872baed6f
# ╠═8771f120-8253-11eb-2e91-d70d9914670e
# ╟─c55c41e0-827e-11eb-0380-d9d33c6e7148
# ╟─d504db00-8311-11eb-2b50-5d6cc93e92a1
# ╠═c26fc67e-8311-11eb-07d3-05a06845d290
# ╟─ad3d8a40-8253-11eb-3709-870c8ba3b76f
# ╠═ad4663e0-8253-11eb-0e02-31d856b9d844
# ╠═f5bb8c76-9068-4065-9f73-edfa80778e95
# ╠═8f953ff7-5d2e-49cf-9abb-697c397c89c6
# ╠═fc1aaaf8-4ec1-485f-ae16-fc21e87212a8
# ╟─ad769b9e-8253-11eb-1946-f7f19dded706
# ╠═ad80d4d0-8253-11eb-36cf-5373bb8741ed
# ╠═91048f97-1725-4266-8adf-18592130cdd2
# ╠═befe81cf-8965-442e-b641-22627770d248
# ╠═9bcf96fe-e5ef-4d26-bdbf-740ac883b82e
# ╟─ada43b50-8253-11eb-2456-87a605b4c2da
# ╠═adc3f852-8253-11eb-3524-c1b3a12dc1e4
# ╠═c65f00f0-7596-483d-8fc0-6f0c4007248f
# ╠═7e83479d-e096-4315-a33e-3299f211b045
# ╠═509f54a9-ef72-40eb-8848-e1467abb5610
# ╟─2d7f9d77-0c90-40b4-87de-74dfe036d7f8
# ╠═2ce19aa2-a158-4724-8c31-9e725096266c
# ╟─cac0286b-1535-436a-bab5-f779cae32998
# ╠═8d017053-31e6-4025-8047-d704e394db08
# ╟─924ccf0c-eb30-4063-8856-d8f36b275ad9
# ╠═1d37024b-3808-4911-9135-5fa53fa4133a
# ╟─f1669dee-8254-11eb-3788-393adbf558ec
# ╟─d9e1dbc2-827e-11eb-160f-038db46faaa6
# ╠═add8b8d0-8253-11eb-0cdc-8f1e28af9e2c
# ╠═b79f2b23-1340-4e58-8adb-d3dd65fbc2c7
# ╠═adf2a972-8253-11eb-05f7-0bc93fd90563
# ╠═66adab50-82f8-11eb-0603-f3a50a03c3b4
# ╠═ab566bc0-82f8-11eb-2b30-b9d62de30027
# ╠═5b1ac8a0-8301-11eb-0e9f-7b4b8bb2b0bd
# ╠═fc17a380-82f8-11eb-30b6-555448b1d939
# ╠═2b3361a0-82f8-11eb-2cf1-a5bfeb3b83f5
# ╟─c0d4b620-82fb-11eb-3735-935581c41b80
# ╠═35f40af0-8315-11eb-2821-efd5d319ffcb
# ╠═b04f5e20-8253-11eb-1263-934f17649375
# ╠═af888bb0-8253-11eb-2dd5-5d0739fd1f15
# ╟─e2fe51e0-8312-11eb-38e0-0d996d94ced1
# ╠═e845a360-8312-11eb-1551-b3e2ecf7cbd2
# ╠═fdd49a2c-5c03-4d2a-9910-7aa3f4a88fc4
# ╠═19387ec0-834f-11eb-347f-69385fa70bd3
# ╟─a49f0c90-8255-11eb-03da-e97868878c21
# ╟─f3d5f660-82fb-11eb-156f-b9443de5879c
# ╠═fb8df9c0-82fb-11eb-0c21-1fb0470f0842
# ╟─ffd1bdf0-82fb-11eb-142d-9bc50fafd5ed
# ╠═af01f690-8253-11eb-232e-bb504f4c3535
# ╠═f1e63450-82fb-11eb-305d-9bec0e0cd50a
# ╟─c39fe410-82fb-11eb-25ed-95925ebc31ba
# ╟─bdd5da00-8277-11eb-08d7-0defcf2c7ba0
