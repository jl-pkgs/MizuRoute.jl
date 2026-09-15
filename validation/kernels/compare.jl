using MizuRoute
using Printf

csv = length(ARGS) >= 1 ? ARGS[1] : "validation/kernels/work/kernel_reference.csv"
lines = readlines(csv)
length(lines) > 1 || error("reference CSV is empty: $csv")

net = RiverNetwork([1], [0])
p = ReachParameters(1;
    length=4500.0,
    slope=0.002,
    bottom_width=15.0,
    side_slope=0.5,
    floodplain_slope=1000.0,
    bankfull_depth=1.0e6,
    mann_n=0.035,
    irf_velocity=1.0,
    irf_diffusivity=1000.0,
)
dt = 3600.0
kw = RoutingModel(EulerKinematicWave(cells_per_reach=20), net, p; dt=dt)
mc = RoutingModel(MuskingumCunge(), net, p; dt=dt)
dw = RoutingModel(DiffusiveWave(cells_per_reach=20), net, p; dt=dt)

maxerr = Dict(
    "kw_qout" => 0.0, "kw_vol" => 0.0,
    "mc_qout" => 0.0, "mc_vol" => 0.0,
    "dw_qout" => 0.0, "dw_vol" => 0.0,
)
worst = Dict(k => 0 for k in keys(maxerr))

function check!(name, got, ref, t)
    err = abs(got - ref)
    if err > maxerr[name]
        maxerr[name] = err
        worst[name] = t
    end
    tol = 5e-9 + 5e-11 * abs(ref)
    err <= tol || error("$name mismatch at step $t: Julia=$got Fortran=$ref abs_error=$err tol=$tol")
end

for line in lines[2:end]
    isempty(strip(line)) && continue
    x = split(line, ',')
    length(x) == 9 || error("unexpected reference row: $line")
    t = parse(Int, strip(x[1]))
    qin = parse(Float64, strip(x[2]))
    qlat = parse(Float64, strip(x[3]))
    ref_kw_q = parse(Float64, strip(x[4])); ref_kw_v = parse(Float64, strip(x[5]))
    ref_mc_q = parse(Float64, strip(x[6])); ref_mc_v = parse(Float64, strip(x[7]))
    ref_dw_q = parse(Float64, strip(x[8])); ref_dw_v = parse(Float64, strip(x[9]))

    got_kw = MizuRoute._route_reach!(kw.method, kw.state, p, 1, qin, qlat, dt)
    got_mc = MizuRoute._route_reach!(mc.method, mc.state, p, 1, qin, qlat, dt)
    got_dw = MizuRoute._route_reach!(dw.method, dw.state, p, 1, qin, qlat, dt)

    check!("kw_qout", got_kw, ref_kw_q, t)
    check!("kw_vol", kw.state.volume[1], ref_kw_v, t)
    check!("mc_qout", got_mc, ref_mc_q, t)
    check!("mc_vol", mc.state.volume[1], ref_mc_v, t)
    check!("dw_qout", got_dw, ref_dw_q, t)
    check!("dw_vol", dw.state.volume[1], ref_dw_v, t)
end

println("mizuRoute Fortran-kernel regression passed")
for name in sort!(collect(keys(maxerr)))
    @printf("%-8s max_abs_error = %.6e at step %d\n", name, maxerr[name], worst[name])
end
