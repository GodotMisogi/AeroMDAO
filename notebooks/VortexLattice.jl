### A Pluto.jl notebook ###
# v0.14.2

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

# ╔═╡ 874269d0-819c-11eb-09c9-1fb09eae02f0
begin
	using Revise
	using AeroMDAO
end

# ╔═╡ 8a9676d0-819c-11eb-1e05-ebf59504b407
using Plots

# ╔═╡ 8c57aed0-819c-11eb-2d1c-0323216e42f0
using StaticArrays

# ╔═╡ 93a864e0-819c-11eb-26d9-71c9b371f684
begin
	using PlutoUI
	alert(text) = Markdown.MD(Markdown.Admonition("warning", "Alert!", [text]))
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));
	definition(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("warning", "Definition", [text]));
	exercise(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("correct", "Exercise", [text]));
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
end;

# ╔═╡ e9306402-5ee0-11eb-1cbb-8dd6f2e08ad6
md"""
# AeroMDAO -- Vortex Lattice Method

**References**: 

1. Mark Drela. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
2. Joseph Katz and Allen Plotkin. _Low-Speed Aerodynamics, Second Edition_. Cambridge University Press, 2001.
"""

# ╔═╡ 80b46780-819c-11eb-1d0a-9fbeb8407ced
md"""
## Your First Wing
"""

# ╔═╡ 80b503be-819c-11eb-18e7-b594665f588b
definition(md"A **wing section** consists of two foil profiles and their chord lengths and twist angles. Between them is their span length with associated _leading-edge_ dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for spans, dihedrals, sweeps, and ``n+1`` entries for foils, chords, and twists for some ``n \in \mathbb N``.")

# ╔═╡ 80b5c710-819c-11eb-38b6-8125b693edc8
md"""![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)"""

# ╔═╡ 80b9e5c0-819c-11eb-0f73-dbc043254a6a
md"""
You can define a wing in this parametrization by using the `HalfWing` type:

```julia
HalfWing(foils  	:: Vector{Foil}, 	# Foil profiles
		 chords 	:: Vector{Real}, 	# Chord lengths
		 twists 	:: Vector{Real}, 	# Twist angles
		 spans  	:: Vector{Real}, 	# Section span lengths
		 dihedrals  :: Vector{Real}, 	# Dihedral angles
		 sweeps 	:: Vector{Real})	# Leading-edge sweep angles
```
"""

# ╔═╡ c0db55c0-830f-11eb-1193-e3275d457f47
foil_path = "..\\data\\airfoil_database\\sd7037.dat";

# ╔═╡ 906fdb90-830f-11eb-1c1c-8f35acf30c59
sd7037 = read_foil(foil_path);

# ╔═╡ ea833240-819c-11eb-299e-216fc8c75464
foils = Foil.([ sd7037, sd7037 ]);

# ╔═╡ 80c5f3b0-819c-11eb-071c-af01954bed9d
wing_right = HalfWing(foils,
                      [1.0, 0.6],
					  [0., 0.],
					  [5.0],
					  [11.3],
					  [2.29]);

# ╔═╡ 80c6de10-819c-11eb-15b5-6fb90139ca04
md"We can create a `Wing` by feeding two `HalfWing`s to it:
```julia
Wing(left 	:: HalfWing, 
	 right 	:: HalfWing)
```
"

# ╔═╡ 80cc8360-819c-11eb-1606-814d961bd0b7
md"In this case, we'd like a symmetric wing, so we just feed `wing_right` as both arguments."

# ╔═╡ 80d0f030-819c-11eb-2ff2-8ff76d3bbb8b
wing = Wing(wing_right, wing_right);

# ╔═╡ 80d1b380-819c-11eb-057a-652abcf73444
md"Now let's see what the outline of our wing looks like, using the following function to get the coordinates:

```julia
plot_wing(some_wing :: Union{Wing, HalfWing})
```
"

# ╔═╡ 80d69580-819c-11eb-0386-d9eab5152a45
begin
	φ_v = @bind φv Slider(0:1e-2:90, default = 15)
	ψ_v = @bind ψv Slider(0:1e-2:90, default = 30)
	z_v = @bind zv Slider(0:1e-2:3, default = 1.)
	md"""
	Horizontal: $(φ_v)
	Vertical: $(ψ_v)
	
	*z*-scale: $(z_v)
	"""
end

# ╔═╡ 80d758d0-819c-11eb-271c-65d473efee7a
gr() # Change to plotlyjs() for interactive animation, but it's slow to render.

# ╔═╡ 80deabce-819c-11eb-1c8b-95f6c533dc06
begin
	plt1 = plot(
				xaxis = "x", yaxis = "y", zaxis = "z",
		 		aspect_ratio = 1, 
		 		camera = (φv, ψv),
		 		zlim = (-0.1, zv)
				)
	
	wing_outline = plot_wing(wing)
	
	plot!(wing_outline, label = "Wing")
end

# ╔═╡ 80e22e3e-819c-11eb-022e-97bf4d5b187a
wing_span, wing_area, wing_mac, wing_AR = span(wing), projected_area(wing), mean_aerodynamic_chord(wing), aspect_ratio(wing);

# ╔═╡ 80e318a0-819c-11eb-1b92-4b3fa6e70be6
md"""
Parameter | Value
:-------- | -----:
Span ($m$)     | $wing_span
Planform Area ($m^2$) | $wing_area 
Mean Aerodynamic Chord ($m$) | $wing_mac
Aspect Ratio | $wing_AR
"""

# ╔═╡ 80e90c10-819c-11eb-1592-1ba9915cac84
exercise(md"Create an anti-symmetric `Wing` consisting of 4 spanwise sections with different NACA 4-digit foil profiles using whatever dimensions you like, as long as they're physically reasonable.")

# ╔═╡ 80eeb160-819c-11eb-207d-d339e232cbb3
md"## Vortex Lattice Analysis

For a 3D case, we use the vortex lattice method for initial designs, given its quick speed for fast analyses. For an analysis, you require ambient reference conditions. In this case, you need the density $\rho$ and a reference location $r_\text{ref}$ for calculating forces and moments."

# ╔═╡ 80ef74b0-819c-11eb-22b2-65555b613023
begin
	ρ = 1.225
	ref = [0.50, 0., 0.]
end;

# ╔═╡ 80fba9b0-819c-11eb-155a-d98d67761735
begin
	U = 10.0
	α = 1.0
	β = 3.0
	Ω = [0.0, 0.0, 0.0]
end;

# ╔═╡ 80f58f30-819c-11eb-2b0a-db7d11ba23bc
md"""

For the boundary conditions, given by the following, we use a `Freestream` type:

Parameter |  Value  | Description
:-------- | :-----: | ----------:
$U_\infty$|    $U    | Freestream speed
$\alpha$  |    $α    | Angle of attack (deg)
$\beta$   |    $β    | Sideslip angle (deg)
$\Omega$  |    $(Ω[1], Ω[2], Ω[3]) | Rotation vector
"""

# ╔═╡ 80fc6d00-819c-11eb-2380-31605f1b2631
freestream = Freestream(U, α, β, Ω);

# ╔═╡ 6982eab0-81a0-11eb-0172-03a207172b04
alert(md"The vortex lattice method only gives reasonable results for small angles of attack and sideslip.")

# ╔═╡ 810546a0-819c-11eb-0232-efe3ec97fc26
md"Now we run the case with specifications of the number of spanwise and chordwise panels by calling the `solve_case()` function, which has an associated method:
```julia
solve_case(wing 				:: Union{Wing, HalfWing},
		   freestream 			:: Freestream, 
		   ρ, 											# Freestream density
		   r_ref = [0.25, 0, 0]; 						# Reference location
		   span_num 			:: Integer,			 	# Number of spanwise panels
		   chord_num 			:: Integer) 			# Number of chordwise panels

```
"

# ╔═╡ 81063100-819c-11eb-15b2-9d20b017ed47
md"It returns nearfield and farfield coefficients, and other arrays for plotting purposes or further analyses."

# ╔═╡ 21111d8e-830f-11eb-2589-dde883171dc9
typeof(wing) <: Wing

# ╔═╡ 810e6e60-819c-11eb-16b2-2943f30a18a1
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, freestream; 
		   rho_ref = ρ, 
		   r_ref = ref, 
		   area_ref = wing_area, 
		   span_ref = wing_span, 
		   chord_ref = wing_mac, 
		   viscous = true, 
		   x_tr = 0.4, 
		   span_num = 48, 
		   chord_num = 12
		  );

# ╔═╡ 810fa6e0-819c-11eb-258f-fb8e4c74369a
HTML(print_coefficients("Wing", nf_coeffs, ff_coeffs; browser = true))

# ╔═╡ 9f7fb000-8306-11eb-1630-53c2be904d85
size(horseshoe_panels)

# ╔═╡ cb9012fa-3234-4635-90f3-9bea0fa7e572
begin
	CL_nf, CL_ff = nf_coeffs[5], ff_coeffs[5]
	CY_nf, CY_ff = nf_coeffs[4], ff_coeffs[4]
	CDi_nf, CDi_ff = nf_coeffs[3], ff_coeffs[3]
end

# ╔═╡ ff10a514-f249-490c-b8ab-4c5c8f441f98
md"Let's try computing the span efficiency factors, which use a parabolic drag polar approximation:
```math
C_{D_i} = \frac{C_L^2 + C_Y^2}{\pi e AR}
```"

# ╔═╡ a093e82d-4f4b-42a8-87dc-43386b77b50d
e_nf = (CL_nf^2) / (π * aspect_ratio(wing) * CDi_nf)

# ╔═╡ e281d210-8316-11eb-3642-efd8f3029a6c
e_ff = (CL_ff^2) / (π * aspect_ratio(wing) * CDi_ff)

# ╔═╡ 6983d510-81a0-11eb-3a6c-81a91699f8f8
exercise(md"Test your wing for different cases and plot the results of the coefficients.")

# ╔═╡ 811db0a0-819c-11eb-2935-e3cd5d1cf036
md"""### Plotting
AeroMDAO provides convenient functions to get quick plots of wings and their solutions using the vortex lattice method.
"""

# ╔═╡ 811ee91e-819c-11eb-3dad-0708c43b3003
begin
	horseshoe_coords = plot_panels(horseshoe_panels[:])
	camber_coords = plot_panels(camber_panels[:])
	# wing_coords = plot_surface(wing);
end;

# ╔═╡ 8126b150-819c-11eb-058d-d947006aa11d
begin
	φ_s = @bind φ Slider(0:1e-2:90, default = 15)
	ψ_s = @bind ψ Slider(0:1e-2:90, default = 30)
	z_s = @bind z_limit Slider(0:1e-2:3, default = 1.)
	stream = @bind stream CheckBox()
	md"""
	Horizontal: $(φ_s)
	Vertical: $(ψ_s)
	
	*z*-scale: $(z_s)
	Streamlines: $stream
	"""
end

# ╔═╡ 8148df50-819c-11eb-2486-95e681f4b1da
begin
	sec_ys = getindex.(midpoint.(horseshoe_panels[1,:]), 2)
	span_CFs = sum(CFs, dims = 1)
	CDis = getindex.(CFs, 1)
	CYs	 = abs.(getindex.(CFs, 2))
	CLs  = getindex.(CFs, 3);
	CL_loadings = 2sum(Γs, dims = 1)[:] / (freestream.V * wing_span)
	
	plot1 = plot(sec_ys, sum(CDis, dims = 1)[:], label = nothing, ylabel = "CDi")
	plot2 = plot(sec_ys, sum(CYs, dims = 1)[:], label = nothing, ylabel = "CY")
	plot3 = plot(sec_ys, CL_loadings, label = nothing, xlabel = "y", ylabel = "CL Loading")
	plot(plot1, plot2, plot3, layout = (3, 1))
end

# ╔═╡ 104865f0-81a1-11eb-39c5-91c32c702cb8
alert(md"The induced drag coefficients computed using nearfield values may be negative at points in the spanwise distribution due to numerical errors or the panel distribution, depending on the case.")

# ╔═╡ 8132bf40-819c-11eb-2ebd-35a6147c3b4c
md"Seed 1"

# ╔═╡ 813dbbc0-819c-11eb-10ed-cf05e3733707
begin
	num_points = 50
	max_z = 0.1
	y = span(wing) / 2 - 0.01
	seed_1 = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))
end;

# ╔═╡ 813e7f10-819c-11eb-36b7-4b0c78e16cc2
md"Seed 2"

# ╔═╡ 81481c02-819c-11eb-0ad8-97014cf886c6
begin
	span_points = 20
	init        = trailing_chopper(ifelse(β == 0, wing_right, wing), span_points) 
	dx, dy, dz  = 0, 0, 1e-3
	seed_2      = [ init .+ Ref([dx, dy, dz])  ;
					init .+ Ref([dx, dy,-dz])  ]
end;

# ╔═╡ 8131d4e0-819c-11eb-1d24-75120ca1247d
streams = ifelse(stream, plot_streams(freestream, seed_2, horseshoes, Γs, 2, 100), nothing)

# ╔═╡ 8127c2c0-819c-11eb-2d3b-1d5cdbb0bac9
begin
	plot(xaxis = "x", yaxis = "y", zaxis = "z",
		 aspect_ratio = 1, 
		 camera = (φ, ψ),
		 zlim = (-0.1, z_limit)
		)
	# plot!.(horseshoe_coords, color = :orange, label = :none)
	plot!.(camber_coords, color = :gray, label = :none)
	scatter!(tupvector(horseshoe_point.(horseshoe_panels))[:], marker = 1, label = :none)
	if stream
		plot!.(streams, color = :green, label = :none)
	end
	plot!()
end

# ╔═╡ 5fc213dc-a94b-44cf-a33c-9cf53d253e10
md"### Stability Derivatives"

# ╔═╡ 52fff6e6-979a-40a5-8b5f-8b6e05677eee
nf, ff, dvs = 
solve_stability_case(wing, freestream; 
					 rho_ref = ρ, 
					 r_ref = ref, 
					 area_ref = wing_area, 
					 span_ref = wing_span, 
					 chord_ref = wing_mac, 
					 viscous = true, 
					 x_tr = [0.4, 0.3, 0.25], 
					 span_num = 24, 
					 chord_num = 12
		   			 );

# ╔═╡ 4a848535-c761-43e4-a5a9-d0e6921edc08
HTML(print_coefficients("Wing", nf, ff; browser = true))

# ╔═╡ 3232bb1e-db45-477e-bf29-3b7ff94159c2
HTML(print_derivatives("Wing", dvs; browser = true))

# ╔═╡ Cell order:
# ╟─e9306402-5ee0-11eb-1cbb-8dd6f2e08ad6
# ╠═874269d0-819c-11eb-09c9-1fb09eae02f0
# ╠═8a9676d0-819c-11eb-1e05-ebf59504b407
# ╠═8c57aed0-819c-11eb-2d1c-0323216e42f0
# ╟─80b46780-819c-11eb-1d0a-9fbeb8407ced
# ╟─80b503be-819c-11eb-18e7-b594665f588b
# ╟─80b5c710-819c-11eb-38b6-8125b693edc8
# ╟─80b9e5c0-819c-11eb-0f73-dbc043254a6a
# ╠═c0db55c0-830f-11eb-1193-e3275d457f47
# ╠═906fdb90-830f-11eb-1c1c-8f35acf30c59
# ╠═ea833240-819c-11eb-299e-216fc8c75464
# ╠═80c5f3b0-819c-11eb-071c-af01954bed9d
# ╟─80c6de10-819c-11eb-15b5-6fb90139ca04
# ╟─80cc8360-819c-11eb-1606-814d961bd0b7
# ╠═80d0f030-819c-11eb-2ff2-8ff76d3bbb8b
# ╟─80d1b380-819c-11eb-057a-652abcf73444
# ╟─80d69580-819c-11eb-0386-d9eab5152a45
# ╠═80d758d0-819c-11eb-271c-65d473efee7a
# ╠═80deabce-819c-11eb-1c8b-95f6c533dc06
# ╠═80e22e3e-819c-11eb-022e-97bf4d5b187a
# ╟─80e318a0-819c-11eb-1b92-4b3fa6e70be6
# ╟─80e90c10-819c-11eb-1592-1ba9915cac84
# ╟─80eeb160-819c-11eb-207d-d339e232cbb3
# ╠═80ef74b0-819c-11eb-22b2-65555b613023
# ╟─80f58f30-819c-11eb-2b0a-db7d11ba23bc
# ╠═80fba9b0-819c-11eb-155a-d98d67761735
# ╠═80fc6d00-819c-11eb-2380-31605f1b2631
# ╟─6982eab0-81a0-11eb-0172-03a207172b04
# ╟─810546a0-819c-11eb-0232-efe3ec97fc26
# ╟─81063100-819c-11eb-15b2-9d20b017ed47
# ╠═21111d8e-830f-11eb-2589-dde883171dc9
# ╠═810e6e60-819c-11eb-16b2-2943f30a18a1
# ╟─810fa6e0-819c-11eb-258f-fb8e4c74369a
# ╠═9f7fb000-8306-11eb-1630-53c2be904d85
# ╠═cb9012fa-3234-4635-90f3-9bea0fa7e572
# ╟─ff10a514-f249-490c-b8ab-4c5c8f441f98
# ╠═a093e82d-4f4b-42a8-87dc-43386b77b50d
# ╠═e281d210-8316-11eb-3642-efd8f3029a6c
# ╟─6983d510-81a0-11eb-3a6c-81a91699f8f8
# ╟─5fc213dc-a94b-44cf-a33c-9cf53d253e10
# ╠═52fff6e6-979a-40a5-8b5f-8b6e05677eee
# ╠═4a848535-c761-43e4-a5a9-d0e6921edc08
# ╠═3232bb1e-db45-477e-bf29-3b7ff94159c2
# ╟─811db0a0-819c-11eb-2935-e3cd5d1cf036
# ╠═811ee91e-819c-11eb-3dad-0708c43b3003
# ╟─8126b150-819c-11eb-058d-d947006aa11d
# ╠═8127c2c0-819c-11eb-2d3b-1d5cdbb0bac9
# ╠═8131d4e0-819c-11eb-1d24-75120ca1247d
# ╠═8148df50-819c-11eb-2486-95e681f4b1da
# ╟─104865f0-81a1-11eb-39c5-91c32c702cb8
# ╟─8132bf40-819c-11eb-2ebd-35a6147c3b4c
# ╠═813dbbc0-819c-11eb-10ed-cf05e3733707
# ╟─813e7f10-819c-11eb-36b7-4b0c78e16cc2
# ╠═81481c02-819c-11eb-0ad8-97014cf886c6
# ╟─93a864e0-819c-11eb-26d9-71c9b371f684
