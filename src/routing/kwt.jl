const _KWT_ALFA = 5.0 / 3.0
const _KWT_MISSING = -9999.0
const _KWT_VERY_SMALL = floatmin(Float64)

@inline _copy_point(x::KWTPoint{Float64}) =
    KWTPoint{Float64}(x.q, x.tentry, x.texit, x.routed, x.qmod)

@inline function _kwt_has_upstream_flow(model::RoutingModel, i::Int)
    @inbounds for j in model.network.upstream[i]
        model.active[j] && return true
    end
    return false
end

"""Time-step average of irregular wave points, ported from `interp_rch`."""
function _kwt_interp_average(told::Vector{Float64}, qold::Vector{Float64},
    t0::Float64, t1::Float64)
    n = length(told)
    n == length(qold) || throw(ArgumentError("KWT interpolation length mismatch"))
    n >= 2 || throw(ArgumentError("KWT interpolation requires two bracketing points"))
    (told[1] <= t0 && told[end] >= t1) ||
        throw(ArgumentError("KWT interpolation bounds do not bracket the time step"))

    ibeg = 1
    for i in 2:n
        if t0 <= told[i]
            ibeg = i
            break
        end
    end
    iend = 1
    for i in 1:n
        if t1 <= told[i]
            iend = i
            break
        end
    end

    # Both target bounds lie between the same two irregular points.
    if t1 < told[ibeg]
        slope = (qold[ibeg] - qold[ibeg - 1]) / (told[ibeg] - told[ibeg - 1])
        q0 = slope * (t0 - told[ibeg - 1]) + qold[ibeg - 1]
        q1 = slope * (t1 - told[ibeg - 1]) + qold[ibeg - 1]
        return 0.5 * (q0 + q1)
    end

    areab = 0.0
    areae = 0.0
    aream = 0.0
    if t0 < told[ibeg]
        slope = (qold[ibeg] - qold[ibeg - 1]) / (told[ibeg] - told[ibeg - 1])
        q0 = slope * (t0 - told[ibeg - 1]) + qold[ibeg - 1]
        areab = (told[ibeg] - t0) * 0.5 * (q0 + qold[ibeg])
    end
    if t1 < told[iend]
        slope = (qold[iend] - qold[iend - 1]) / (told[iend] - told[iend - 1])
        q1 = slope * (t1 - told[iend - 1]) + qold[iend - 1]
        areae = (t1 - told[iend - 1]) * 0.5 * (qold[iend - 1] + q1)
    end
    if ibeg < iend
        for i in (ibeg + 1):iend
            if i < iend || (i == iend && t1 == told[iend] && t0 < told[iend - 1])
                aream += (told[i] - told[i - 1]) * 0.5 * (qold[i - 1] + qold[i])
            end
        end
    end
    return (areab + areae + aream) / (t1 - t0)
end

"""Port of mizuRoute `remove_rch`: remove the least informative wave points."""
function _kwt_reduce(q::Vector{Float64}, ti::Vector{Float64}, tx::Vector{Float64}, maxq::Int)
    length(q) == length(ti) == length(tx) || throw(ArgumentError("KWT reduction length mismatch"))
    keep = collect(eachindex(q))
    # `q` includes the zero/last-routed point. Fortran exits when MPRT<MAXQPAR,
    # where MPRT excludes that zero point.
    while length(keep) - 1 >= maxq
        bestpos = 0
        besterr = Inf
        for pos in 2:(length(keep) - 1)
            a, b, c = keep[pos - 1], keep[pos], keep[pos + 1]
            den = ti[c] - ti[a]
            den == 0.0 && continue
            qint = q[a] + (q[c] - q[a]) / den * (ti[b] - ti[a])
            err = abs(qint - q[b])
            if err < besterr
                besterr = err
                bestpos = pos
            end
        end
        bestpos == 0 && break
        deleteat!(keep, bestpos)
    end
    return q[keep], ti[keep], tx[keep]
end

"""Port of mizuRoute `kinwav_rch` including kinematic-shock merging."""
function _kwt_kinwav(qin::Vector{Float64}, tin::Vector{Float64},
    reach_length::Float64, slope::Float64, mann_n::Float64,
    tstart::Float64, tend::Float64)
    ni = length(qin)
    ni == length(tin) || throw(ArgumentError("KWT wave/time length mismatch"))
    ni == 0 && return Float64[], Float64[], Float64[], Bool[]

    kcoef = sqrt(slope) / mann_n
    q0 = copy(qin); q1 = copy(qin); q2 = copy(qin)
    t0 = copy(tin); t1 = copy(tin)
    mf = collect(1:ni)
    ix = collect(1:ni)
    wc = zeros(Float64, ni)
    @inbounds for i in 1:ni
        wc[i] = _KWT_ALFA * kcoef^(1.0 / _KWT_ALFA) *
                q1[i]^((_KWT_ALFA - 1.0) / _KWT_ALFA)
    end

    nn = ni
    if nn > 1
        x = 0.0
        while true
            xb = reach_length
            ixb = 0
            for iw in 2:nn
                jw = iw - 1
                (wc[iw] == 0.0 || wc[jw] == 0.0) && continue
                wdiff = 1.0 / wc[jw] - 1.0 / wc[iw]
                wdiff == 0.0 && continue
                xxb = (t1[iw] - t1[jw]) / wdiff
                (xxb < x || xxb > xb) && continue
                xb = xxb
                ixb = iw
            end
            xb == reach_length && break
            ixb == 0 && break

            nn -= 1
            jxb = ixb - 1
            q2[jxb] = max(q2[jxb], q2[ixb])
            q1[jxb] = min(q1[jxb], q1[ixb])
            a2 = (q2[jxb] / kcoef)^(1.0 / _KWT_ALFA)
            a1 = (q1[jxb] / kcoef)^(1.0 / _KWT_ALFA)
            cm = (q2[jxb] - q1[jxb]) / (a2 - a1)
            t1[jxb] = t1[jxb] + xb / wc[jxb] - xb / cm
            wc[jxb] = cm

            first_original = ix[ixb]
            @inbounds for j in first_original:ni
                mf[j] -= 1
            end
            @inbounds for j in ixb:nn
                ix[j] = ix[j + 1]
                t1[j] = t1[j + 1]
                wc[j] = wc[j + 1]
                q1[j] = q1[j + 1]
                q2[j] = q2[j + 1]
            end
            x = xb
        end
    end

    qout = Float64[]
    tout = Float64[]
    texit = Float64[]
    routed = Bool[]
    function addpoint(qnew::Float64, told::Float64, tnew::Float64)
        if !isempty(texit) && tnew <= texit[end]
            tnew = texit[end] + 1.0
        elseif isempty(texit) && tnew <= tstart
            tnew = tstart + 1.0
        end
        push!(qout, qnew)
        push!(tout, told)
        push!(texit, tnew)
        push!(routed, tnew < tend)
    end

    for iroute in 1:nn
        wc[iroute] < _KWT_VERY_SMALL &&
            throw(ArgumentError("zero KWT wave celerity"))
        tx = reach_length / wc[iroute] + t1[iroute]
        tnext = iroute < nn ? reach_length / wc[iroute + 1] + t1[iroute + 1] : floatmax(Float64)
        if q1[iroute] != q2[iroute]
            if tx < tend
                tx2 = min(tx + 1.0, tx + 0.5 * (min(tnext, tend) - tx))
                tx2 == tx && throw(ArgumentError("KWT shock exit times are identical"))
                addpoint(q1[iroute], t1[iroute], tx)
                addpoint(q2[iroute], t1[iroute], tx2)
            else
                for j in 1:ni
                    mf[j] == iroute && addpoint(q0[j], t0[j], tx)
                end
            end
        else
            addpoint(q1[iroute], t1[iroute], tx)
        end
    end
    return qout, tout, texit, routed
end

"""Merge basin and routed-wave series from all immediate upstream reaches.

This is a direct structural port of `qexmul_rch`: basin discharge contributes a
two-point series at T0/T1, routed reach waves contribute their routed points,
and all series are interpolated onto the union of wave exit times after scaling
by upstream/downstream channel width.
"""
function _kwt_qexmul!(model::RoutingModel, i::Int, qlat::Vector{Float64},
    t0::Float64, t1::Float64)
    ups = model.network.upstream[i]
    isempty(ups) && return Float64[], Float64[]
    dt = t1 - t0
    nupb = length(ups)
    nupr = count(j -> _kwt_has_upstream_flow(model, j), ups)
    nups = nupb + nupr

    if nups == 1
        j = ups[1]
        return [qlat[j] / model.params.bottom_width[i]], [t1]
    end

    series = Vector{Vector{KWTPoint{Float64}}}()
    widths = Float64[]

    # Basin routed flow from every immediate upstream segment.
    for j in ups
        push!(series, KWTPoint{Float64}[
            KWTPoint(model.state.qlat_prev[j], t0, t0, true, _KWT_MISSING),
            KWTPoint(qlat[j], t1, t1, true, _KWT_MISSING),
        ])
        push!(widths, 1.0)
    end

    # Reach-wave contribution from non-headwater upstream segments.
    for j in ups
        _kwt_has_upstream_flow(model, j) || continue
        wave = model.state.waves[j]
        isempty(wave) && continue
        nr = count(x -> x.routed, wave)
        ns = length(wave)
        nq = min(nr + 1, ns)
        push!(series, [_copy_point(wave[k]) for k in 1:nq])
        push!(widths, model.params.bottom_width[j])

        # Fortran retains the last routed point plus all non-routed points.
        firstkeep = max(nr, 1)
        model.state.waves[j] = [_copy_point(wave[k]) for k in firstkeep:ns]
    end

    nseries = length(series)
    nseries == 0 && return Float64[], Float64[]
    itim = fill(2, nseries)       # Fortran index 1 (after zero point).
    done = falses(nseries)
    ctime = [s[2].texit for s in series]
    qd = Float64[]
    td = Float64[]
    time_old = -floatmax(Float64)

    while count(done) < nseries
        jups = argmin(ctime)
        if done[jups]
            ctime[jups] = floatmax(Float64)
            continue
        end
        smin = series[jups]
        iw = itim[jups]
        if !smin[iw].routed
            done[jups] = true
            ctime[jups] = floatmax(Float64)
            continue
        end

        target = ctime[jups]
        if target != time_old
            qagg = 0.0
            for k in 1:nseries
                s = series[k]
                ik = itim[k]
                scfac = widths[k] / model.params.bottom_width[i]
                if k == jups
                    flow = s[ik].q
                else
                    ibeg = ik
                    if s[ibeg].texit >= target
                        ibeg -= 1
                    end
                    iend = ibeg + 1
                    (ibeg < 1 || iend > length(s)) &&
                        throw(ArgumentError("KWT upstream interpolation is not bracketed"))
                    ta, tb = s[ibeg].texit, s[iend].texit
                    (tb < target || ta > target) &&
                        throw(ArgumentError("KWT upstream wave times are not ordered"))
                    flow = s[ibeg].q + (s[iend].q - s[ibeg].q) /
                           (tb - ta) * (target - ta)
                end
                qagg += flow * scfac
            end
            push!(qd, qagg)
            push!(td, target)
            time_old = target
        end

        if itim[jups] == length(series[jups])
            done[jups] = true
            ctime[jups] = floatmax(Float64)
        else
            itim[jups] += 1
            ctime[jups] = series[jups][itim[jups]].texit
        end
    end
    return qd, td
end

function _kwt_getusq!(model::RoutingModel, i::Int, qlat::Vector{Float64},
    t0::Float64, t1::Float64)
    qd, td = _kwt_qexmul!(model, i, qlat, t0, t1)
    isempty(qd) && throw(ArgumentError("KWT received no upstream wave points"))
    dt = t1 - t0
    wave = model.state.waves[i]
    if isempty(wave)
        wave = KWTPoint{Float64}[
            KWTPoint(qd[1], t0 - dt, t0, true, _KWT_MISSING)
        ]
    else
        wave = [_copy_point(x) for x in wave]
    end

    q = [x.q for x in wave]
    ti = [x.tentry for x in wave]
    tx = [x.texit for x in wave]
    append!(q, qd)
    append!(ti, td)
    append!(tx, fill(_KWT_MISSING, length(qd)))
    return q, ti, tx
end

"""Exact mizuRoute-style KWT network step.

KWT cannot be expressed as the generic one-reach flux update: routed wave points
are removed from upstream state while the downstream reach is processed, which
is what permits one wave to cross several reaches within a single model step.
"""
function step!(model::RoutingModel{M,S,T}, qlat_in::AbstractVector;
    water_management=nothing) where {M<:LagrangianKWT,S<:KWTState,T}
    n = length(model.network)
    length(qlat_in) == n || throw(ArgumentError("qlat length mismatch"))
    if water_management !== nothing
        length(water_management) == n || throw(ArgumentError("water_management length mismatch"))
        any(x -> Float64(x) != 0.0, water_management) &&
            throw(ArgumentError("exact KWT water-management wave modification is not implemented yet"))
    end

    qlat = max.(Float64.(qlat_in), 0.0)
    t0 = model.time
    t1 = t0 + model.dt
    state = model.state

    for i in model.network.order
        if !_kwt_has_upstream_flow(model, i)
            # Official KWT does not route within a headwater reach; basin-routed
            # runoff is emitted directly and the wave structure is reset.
            state.qin[i] = 0.0
            state.qout[i] = qlat[i]
            state.volume[i] = 0.0
            state.wm_actual[i] = 0.0
            state.waves[i] = KWTPoint{Float64}[
                KWTPoint(_KWT_MISSING, _KWT_MISSING, _KWT_MISSING, false, _KWT_MISSING)
            ]
            continue
        end

        q_up = 0.0
        for j in model.network.upstream[i]
            model.active[j] && (q_up += state.qout[j])
        end
        state.qin[i] = q_up
        state.wm_actual[i] = 0.0
        state.volume[i] = 0.0

        q, ti, tx = _kwt_getusq!(model, i, qlat, t0, t1)
        if length(q) > model.method.max_packets
            q, ti, tx = _kwt_reduce(q, ti, tx, model.method.max_packets)
        end

        # Exclude zero/last-routed point for kinwav_rch.
        qr, tir, txr, routed = _kwt_kinwav(q[2:end], ti[2:end],
            model.params.length[i], model.params.slope[i], model.params.mann_n[i], t0, t1)
        qj = vcat(q[1], qr)
        tij = vcat(ti[1], tir)
        txj = vcat(tx[1], txr)
        fr = vcat(true, routed)
        nr = count(fr) - 1
        nq2 = length(qr)
        nr + 2 <= length(qj) || throw(ArgumentError("KWT lacks next wave for time averaging"))

        qavg = _kwt_interp_average(txj[1:(nr + 2)], qj[1:(nr + 2)], t0, t1)
        state.qout[i] = qavg * model.params.bottom_width[i] + qlat[i]

        # Insert the interpolated point at T_END exactly as Fortran does.
        ia = nr + 1
        ib = nr + 2
        den = txj[ib] - txj[ia]
        den == 0.0 && throw(ArgumentError("KWT interpolation has duplicate exit times"))
        frac = (t1 - txj[ia]) / den
        qend = qj[ia] + (qj[ib] - qj[ia]) * frac
        timei = tij[ia] + (tij[ib] - tij[ia]) * frac

        newwave = KWTPoint{Float64}[]
        for k in 1:(nr + 1)
            push!(newwave, KWTPoint(qj[k], tij[k], txj[k], fr[k], _KWT_MISSING))
        end
        push!(newwave, KWTPoint(qend, timei, t1, true, _KWT_MISSING))
        for k in (nr + 2):(nq2 + 1)
            push!(newwave, KWTPoint(qj[k], tij[k], txj[k], fr[k], _KWT_MISSING))
        end
        state.waves[i] = newwave

        # At an outlet, routed points are immediately discarded, leaving the
        # interpolated T_END point as the new zero element.
        if model.network.downstream[i] == 0
            firstkeep = nr + 2
            state.waves[i] = [_copy_point(newwave[k]) for k in firstkeep:length(newwave)]
        end
    end

    state.qlat_prev .= qlat
    model.time = t1
    return state.qout
end

# The exact KWT path is network-coupled and therefore intentionally has no
# `_route_reach!` implementation. Calling it directly would omit upstream wave
# extraction/merging and is a model error.
function _route_reach!(::LagrangianKWT, ::KWTState, ::ReachParameters,
    ::Int, ::Float64, ::Float64, ::Float64)
    throw(ArgumentError("LagrangianKWT must be advanced with step!, not a reach-local kernel"))
end
