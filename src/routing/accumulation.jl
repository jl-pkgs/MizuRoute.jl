function _route_reach!(::Accumulation, state::BasicState, p::ReachParameters,
    i::Int, qin::Float64, qlat::Float64, dt::Float64)
    qout = qin + qlat
    state.qout[i] = qout
    state.volume[i] = 0.0
    return qout
end
