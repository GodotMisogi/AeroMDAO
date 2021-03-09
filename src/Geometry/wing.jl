#-----------------WING---------------------#

"""
	HalfWing(foils      :: Vector{Foil},
			 chords     :: Vector{Real},
			 twists     :: Vector{Real},
			 spans      :: Vector{Real},
			 dihedrals  :: Vector{Real},
			 sweeps     :: Vector{Real})
			 
Definition for a `HalfWing` consisting of ``N+1`` airfoils and their associated chord lengths ``c``, twist angles ``\\iota``, for ``N`` sections with span lengths ``b``, dihedrals ``\\delta`` and sweep angles ``\\Lambda``, with all angles in degrees.
"""
struct HalfWing{T <: Real} <: Aircraft
	foils       :: Vector{Foil{T}}
	chords      :: Vector{T}
	twists      :: Vector{T}
	spans       :: Vector{T}
	dihedrals   :: Vector{T}
	sweeps      :: Vector{T}
	HalfWing(foils :: Vector{Foil{T}}, chords :: Vector{T}, twists :: Vector{T}, spans :: Vector{T}, dihedrals :: Vector{T}, sweeps :: Vector{T}) where T <: Real = new{T}(foils, chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps))
end

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
quarter_chord(chord) = 0.25 * chord

"""
	span(half_wing :: HalfWing)

Compute the planform span of a `HalfWing`.
"""
span(wing :: HalfWing) = sum(wing.spans)

"""
	projected_area(half_wing :: HalfWing)
	
Compute the projected area of a `HalfWing`.
"""
function projected_area(wing :: HalfWing)
	mean_chords = fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
	mean_twists = fwdsum(wing.twists) / 2 # Mean twist angles of sections
	sum(@. wing.spans * mean_chords * cos(mean_twists))
end

"""
	mean_aerodynamic_chord(half_wing :: HalfWing)
	
Compute the mean aerodynamic chord of a `HalfWing`.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
	mean_chords = fwdsum(wing.chords) / 2
	taper_ratios = fwddiv(wing.chords)
	areas = mean_chords .* wing.spans
	macs = mean_aerodynamic_chord.(wing.chords[1:end-1], taper_ratios)
	sum(macs .* areas) / sum(areas)
end

"""
	chop_sections(set1, set2, n :: Integer; spacing = "cosine")

Divide two vectors into ``n`` sections, with a named argument for spacing.
"""
function chop_sections(set1, set2, n :: Integer; spacing = "uniform") 
	if spacing == "cosine"
		space = cosine_dist(0.5, 1., n + 1)
	else
		space = range(0, stop = 1, length = n + 1)
	end

	[ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

coords_chopper(coords, n, spacing = "uniform") = [ chop_sections.(coords[1:end-1], coords[2:end], n; spacing = spacing)...; [coords[end]] ]
chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) ∈ zip(lead, trail) ]
chord_chopper(coords, n) = [ [ weighted_vector(chord[1,:], chord[2,:], μ) for μ ∈ cosine_dist(0.5, 1., n + 1) ] for chord ∈ coords ]
span_chopper(lead, trail, div) = coords_chopper(lead, div), coords_chopper(trail, div)
wing_chopper(lead, trail, span_num, chord_num) = chord_chopper(chord_sections(span_chopper(lead, trail, span_num)...), chord_num)

"""
	leading_edge(wing :: HalfWing, flip = false)

Compute the leading edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function leading_edge(wing :: HalfWing, flip = false)
	spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

	sweeped_spans = [ 0; cumsum(@. spans * tan(sweeps)) ]
	dihedraled_spans = [ 0; cumsum(@. spans * tan(dihedrals)) ]
	cum_spans = [ 0; cumsum(spans) ]
	
	leading = SVector.(sweeped_spans, ifelse(flip, -cum_spans, cum_spans), dihedraled_spans)

	ifelse(flip, leading[end:-1:1], leading)
end

leading_chopper(obj :: HalfWing, span_num; flip = false) = coords_chopper(leading_edge(obj, flip), span_num)
trailing_chopper(obj :: HalfWing, span_num; flip = false) = coords_chopper(leading_edge(obj, flip), span_num)

"""
	wing_bounds(wing :: HalfWing, flip = false)

Compute the leading and trailing edge coordinates of a `HalfWing`, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip = false)
	chords = wing.chords
	twisted_chords = @. chords * sin(wing.twists)
	
	leading = leading_edge(wing, flip)
	trailing = SVector.(chords, (zeros ∘ length)(chords), twisted_chords) 
	shifted_trailing = ifelse(flip, trailing[end:-1:1], trailing) .+ leading

	leading, shifted_trailing
end

"""
	wing_coords(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the coordinates of a `HalfWing` consisting of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function wing_coords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer, flip = false)
	leading_xyz = leading_edge(wing, flip)
	temp_foils = @. wing.chords * (extend_yz ∘ cosine_foil)(wing.foils, chord_num)
	
	scaled_foils = ifelse(flip, temp_foils[end:-1:1], temp_foils)
	twists = ifelse(flip, wing.twists[end:-1:1], wing.twists)

	foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

	coords_chopper(foil_coords, span_num)
end

"""
	camber_coords(wing :: HalfWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the coordinates of a `HalfWing` consisting of camber distributions of `Foil`s and relevant geometric quantities, given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
function camber_coords(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer, flip = false)
	leading_xyz = leading_edge(wing, flip)
	
	scaled_foils = @. wing.chords * (camber_coordinates ∘ camber_thickness)(wing.foils, chord_num)

	scaled_foils = ifelse(flip, scaled_foils[end:-1:1], scaled_foils)
	twists = ifelse(flip, wing.twists[end:-1:1], wing.twists)

	foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

	coords_chopper(foil_coords, span_num)
end

"""
	make_panels(coords)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
function make_panels(coords) 
	spanlist = vectarray.(coords)    
	secs1secs2 = zip(spanlist, spanlist[2:end])
	hcat(( Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) ∈ secs1secs2 )...)
end


"""
	mesh_horseshoes(obj :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for lifting-line/vortex lattice analyses using horseshoe elements.
"""
mesh_horseshoes(obj :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; flip = false) = (make_panels ∘ wing_chopper)(wing_bounds(obj, flip)..., span_num, chord_num)

"""
	mesh_wing(wing :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for 3D analyses using doublet-source elements or equivalent formulations. TODO: Tip meshing.
"""
mesh_wing(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; flip = false) = (make_panels ∘ wing_coords)(wing, span_num, chord_num, flip)

"""
	mesh_cambers(wing :: HalfWing, n_s :: Integer, n_c :: Integer; flip = false)

Mesh the camber distribution of a `HalfWing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for vortex lattice analyses.
"""
mesh_cambers(wing :: HalfWing, span_num :: Vector{<: Integer}, chord_num :: Integer; flip = false) = (make_panels ∘ camber_coords)(wing, span_num, chord_num, flip)

"""
	Wing(left :: HalfWing, right :: HalfWing)

A composite type consisting of two `HalfWing`s with fields `left` and `right` for constructing a wing.
"""
struct Wing{T <: Real} <: Aircraft
	left :: HalfWing{T}
	right :: HalfWing{T}
end

"""
	span(wing :: Wing)
	
Compute the span of a `Wing`.
"""
span(wing :: Wing) = span(wing.left) + span(wing.right)

"""
	projected_area(wing :: Wing)
	
Compute the projected_area of a `Wing`.
"""
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)

"""
	mean_aerodynamic_chord(wing :: Wing)
	
Compute the mean aerodynamic chord of a `Wing`.
"""
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2

"""
	mean_aerodynamic_chord(wing :: Union{Wing, HalfWing})
	
Compute the mean aerodynamic chord of a `HalfWing` or `Wing`.
"""
aspect_ratio(wing :: Union{Wing, HalfWing}) = aspect_ratio(span(wing), projected_area(wing))

"""
	wing_bounds(wing :: Wing)

Return the leading and trailing edge coordinates of a `Wing`.
"""
function wing_bounds(wing :: Wing)
	left_lead, left_trail 	= wing_bounds(wing.left, true)
	right_lead, right_trail = wing_bounds(wing.right)

	leading 	= [ left_lead; right_lead ]
	trailing 	= [ left_trail; right_trail ]

	leading, trailing
end

leading_chopper(obj :: Wing, span_num :: Integer; flip = false) = span_chopper(wing_bounds(obj)..., span_num)[1]

trailing_chopper(obj :: Wing, span_num :: Integer; flip = false) = span_chopper(wing_bounds(obj)..., span_num)[2]

"""
	mesh_horseshoes(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for lifting-line/vortex lattice analyses using horseshoe elements.
"""
function mesh_horseshoes(wing :: Wing, span_num, chord_num)
	left_panels = mesh_horseshoes(wing.left, span_num[end:-1:1], chord_num, flip = true)
	right_panels = mesh_horseshoes(wing.right, span_num, chord_num)

	[ left_panels right_panels ] 
end

"""
	mesh_wing(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for 3D analyses using doublet-source elements or equivalent formulations. TODO: Tip meshing.
"""
function mesh_wing(wing :: Wing, span_num, chord_num)
	left_panels = mesh_wing(wing.left, span_num[end:-1:1], chord_num, flip = true)
	right_panels = mesh_wing(wing.right, span_num, chord_num)

	[ left_panels right_panels ] 
end

"""
	mesh_cambers(wing :: Wing, n_s :: Integer, n_c :: Integer)

Mesh the camber distribution of a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions meant for vortex lattice analyses.
"""
function mesh_cambers(wing :: Wing, span_num, chord_num)
	left_panels = mesh_cambers(wing.left, span_num[end:-1:1], chord_num, flip = true)
	right_panels = mesh_cambers(wing.right, span_num, chord_num)

	[ left_panels right_panels ] 
end

"""
	vlmesh_wing(wing :: Union{Wing, HalfWing}, n_s :: Integer, n_c :: Integer)

Return the horseshoe and camber panel distributions of a `HalfWing` or `Wing` consisting of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions for vortex lattice analyses.
"""
vlmesh_wing(wing :: Union{Wing, HalfWing}, span_num, chord_num) = mesh_horseshoes(wing, span_num, chord_num), mesh_cambers(wing, span_num, chord_num)

function paneller(wing :: Union{Wing, HalfWing}, span_num, chord_num; rotation = one(RotMatrix{3, Float64}), translation = SVector(0,0,0))
	horseshoes, cambers = vlmesh_wing(wing, span_num, chord_num) 
	[ transform(panel, rotation, translation) for panel in horseshoes ],
	[ transform(panel, rotation, translation) for panel in cambers ]
end

"""
	info(wing :: Union{Wing, HalfWing})

Return the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
info(wing :: Union{Wing, HalfWing}) = span(wing), projected_area(wing), mean_aerodynamic_chord(wing), aspect_ratio(wing)

"""
	print_info(wing :: Union{Wing, HalfWing})

Print the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
function print_info(wing :: Union{Wing, HalfWing})
	span, area, mac, AR = info(wing)
	println("Span: ", span, " m")
	println("Area: ", area, " m²")
	println("MAC: ", mac, " m")
	println("Aspect Ratio: ", AR)
end
