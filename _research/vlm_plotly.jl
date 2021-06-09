##
using PlotlyJS

# layout = Layout(
#                 title = "Vortex Lattice Aircraft",
#                 scene=attr(aspectmode="manual", aspectratio=attr(x=1,y=1,z=1)),
#                 # zlim=(-0.1, 5.0)
#                 )

span_points = 20
wing_init   = trailing_chopper(wing.right, span_points) 
htail_init  = leading_chopper(htail.right, span_points)

dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz])  ; 
                init .+ Ref([dx, dy, -dz]) ;
                htail_init .+ Ref(htail_location .+ [dx, dy, dz])];


layout, plt = plot_case(horseshoe_panels, camber_panels, Γs, horseshoes, freestream, seed, 2, 100);

##
plt

##
savefig(plt, "VLM.pdf", format = "pdf", width = 1920, height = 1080, scale = 5)

## Streamlines
# reset_timer!()

# @timeit "Computing Streamlines" trace_streams = trace_streamlines(freestream, horseshoe_panels, horseshoes, Γs, 2, 100);

# print_timer()

# trace_horses = trace_panels(horseshoe_panels, Γs)
# trace_cambers = trace_panels(camber_panels)

# plot([ 
#         [ trace for trace in trace_cambers ]...,
#         # [ trace for trace in trace_streams ]...,
#         [ trace for trace in trace_horses ]...,
#      ], 
#      layout)