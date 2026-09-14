function _check_model_inputs(net::RiverNetwork, p::ReachParameters, dt::Real)
    length(net) == length(p.length) || throw(ArgumentError("network and ReachParameters have different reach counts"))
    dt > 0 || throw(ArgumentError("routing dt must be positive"))
    return nothing
end

function _initial_state(::Accumulation, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    BasicState(fill(q0, n), zeros(n), zeros(n), zeros(n))
end

function _initial_state(method::IRF, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    kernels = Vector{Vector{Float64}}(undef, n)
    history = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        kernels[i] = _irf_weights(p.length[i], p.irf_velocity[i], p.irf_diffusivity[i], dt;
            horizon_factor=method.horizon_factor, min_steps=method.min_steps)
        history[i] = fill(q0, length(kernels[i]))
    end
    IRFState(fill(q0, n), zeros(n), zeros(n), zeros(n), history, kernels)
end

function _initial_state(::LagrangianKWT, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    packets = [WavePacket{Float64}[] for _ in 1:n]
    KWTState(fill(q0, n), zeros(n), zeros(n), zeros(n), packets)
end

function _initial_state(method::EulerKinematicWave, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    method.cells_per_reach >= 3 || throw(ArgumentError("Euler KW requires cells_per_reach >= 3"))
    profiles = [fill(q0, method.cells_per_reach + 1) for _ in 1:n]
    EulerKWState(fill(q0, n), zeros(n), zeros(n), zeros(n), profiles)
end

function _initial_state(::MuskingumCunge, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    MCState(fill(q0, n), zeros(n), zeros(n), zeros(n), fill(q0, n), fill(q0, n))
end

function _initial_state(method::DiffusiveWave, net::RiverNetwork, p::ReachParameters, dt::Float64, q0::Float64)
    n = length(net)
    method.cells_per_reach >= 3 || throw(ArgumentError("DiffusiveWave requires cells_per_reach >= 3"))
    0.0 <= method.alpha <= 1.0 || throw(ArgumentError("DiffusiveWave alpha must be in [0,1]"))
    0.0 <= method.beta <= 1.0 || throw(ArgumentError("DiffusiveWave beta must be in [0,1]"))
    profiles = [fill(q0, method.cells_per_reach + 1) for _ in 1:n]
    DWState(fill(q0, n), zeros(n), zeros(n), zeros(n), profiles)
end

"""Construct a routing model with a fixed routing interval `dt` [s]."""
function RoutingModel(method::AbstractRoutingMethod, net::RiverNetwork, p::ReachParameters;
    dt::Real, initial_discharge::Real=0.0)
    _check_model_inputs(net, p, dt)
    Δt = Float64(dt)
    q0 = max(Float64(initial_discharge), 0.0)
    state = _initial_state(method, net, p, Δt, q0)
    RoutingModel(method, net, p, Δt, 0.0, state)
end

@inline discharge(model::RoutingModel) = model.state.qout
@inline inflow(model::RoutingModel) = model.state.qin
@inline reach_storage(model::RoutingModel) = model.state.volume

function _upstream_flow(model::RoutingModel, i::Int)
    q = 0.0
    @inbounds for j in model.network.upstream[i]
        q += model.state.qout[j]
    end
    return q
end

"""Advance all reaches by one routing interval.

Arguments
---------
- `qlat[i]`: local lateral inflow to reach `i` [m³/s].
- `water_management[i]` (optional): positive abstraction and negative injection
  [m³/s], following the mizuRoute sign convention.

The network is always traversed from headwaters to outlets, so downstream
reaches receive the just-computed upstream discharge for the current step.
"""
function step!(model::RoutingModel, qlat::AbstractVector; water_management=nothing)
    n = length(model.network)
    length(qlat) == n || throw(ArgumentError("qlat length $(length(qlat)) != $n reaches"))
    water_management === nothing || length(water_management) == n ||
        throw(ArgumentError("water_management length mismatch"))
    dt = model.dt

    for i in model.network.order
        q_up_raw = _upstream_flow(model, i)
        q_local_raw = max(Float64(qlat[i]), 0.0)
        request = water_management === nothing ? 0.0 : Float64(water_management[i])
        q_up, q_local, v_after_wm, actual = _apply_water_management(
            q_up_raw, q_local_raw, model.state.volume[i], request, dt)
        model.state.volume[i] = v_after_wm
        model.state.qin[i] = q_up_raw
        model.state.wm_actual[i] = actual
        _route_reach!(model.method, model.state, model.params, i, q_up, q_local, dt)
    end
    model.time += dt
    return model.state.qout
end

"""Route a complete lateral-inflow series.

`qlat` must have shape `(nreach, ntime)`. The returned matrix has the same
shape and stores routed reach discharge [m³/s].
"""
function route_series(model::RoutingModel, qlat::AbstractMatrix; water_management=nothing)
    nreach, ntime = size(qlat)
    nreach == length(model.network) || throw(ArgumentError("first dimension must equal number of reaches"))
    q = Matrix{Float64}(undef, nreach, ntime)
    if water_management !== nothing
        size(water_management) == size(qlat) || throw(ArgumentError("water_management matrix size mismatch"))
    end
    for t in 1:ntime
        wm = water_management === nothing ? nothing : view(water_management, :, t)
        step!(model, view(qlat, :, t); water_management=wm)
        q[:, t] .= discharge(model)
    end
    return q
end

"""Reset time and dynamic routing states while retaining network, parameters and method."""
function reset!(model::RoutingModel; initial_discharge::Real=0.0)
    q0 = max(Float64(initial_discharge), 0.0)
    model.state = _initial_state(model.method, model.network, model.params, model.dt, q0)
    model.time = 0.0
    return model
end
