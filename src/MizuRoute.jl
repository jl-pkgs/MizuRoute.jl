module MizuRoute

include("types.jl")
include("network.jl")
include("hydraulics.jl")
include("remap.jl")
include("hillslope.jl")
include("balance.jl")
include("routing/accumulation.jl")
include("routing/irf.jl")
include("routing/kwt.jl")
include("routing/advection_diffusion.jl")
include("routing/euler_kw.jl")
include("routing/muskingum_cunge.jl")
include("routing/diffusive_wave.jl")
include("model.jl")

export RiverNetwork, ReachParameters, RoutingModel
export AbstractRoutingMethod, Accumulation, IRF, LagrangianKWT,
       EulerKinematicWave, MuskingumCunge, DiffusiveWave
export top_width, wetted_perimeter, flow_area, water_height,
       hydraulic_radius, manning_discharge, flow_depth, storage,
       celerity, diffusivity
export gamma_uh_weights, GammaUHRouter, route_hillslope!
export runoff_depth_to_volume, map_runoff!, map_runoff
export step!, route_series, discharge, inflow, reach_storage,
       water_balance_error, reset!
export upstream_indices, downstream_index, headwaters, outlets

end
