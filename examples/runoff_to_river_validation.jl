using MizuRoute
using Printf

function main()
    # 1. 输入：8 个网格的径流深率 [m/s]
    dt = 3600.0
    storm_steps = 12
    cell_area = fill(25e6, 8)
    cell_to_reach = [1, 1, 2, 2, 3, 3, 4, 4]

    # 2. 河网与河段参数
    net = RiverNetwork([101, 102, 103, 104], [103, 103, 104, 0])
    p = ReachParameters(length(net);
        length=[4000.0, 3500.0, 8000.0, 12000.0],
        slope=[0.0020, 0.0030, 0.0012, 0.0008],
        bottom_width=[8.0, 6.0, 15.0, 22.0],
        side_slope=1.0,
        floodplain_slope=20.0,
        bankfull_depth=[2.0, 1.8, 2.8, 3.5],
        mann_n=0.035,
    )

    # 3. 坡面汇流：每个网格分别通过 gamma 单位线
    hillslope = GammaUHRouter(length(cell_area), dt; shape=2.5, scale=3600.0)
    ntime = storm_steps + length(hillslope.kernel) - 1

    # 径流深输入
    runoff_depth = zeros(length(cell_area), ntime)
    runoff_depth[:, 1:storm_steps] .= 1e-6

    # 4. 河道汇流：扩散波
    river = RoutingModel(DiffusiveWave(cells_per_reach=12), net, p; dt=dt)
    qinst = zeros(length(net))
    routed_depth = zeros(length(cell_area))
    qlat = similar(qinst)
    outlet_q = zeros(ntime)
    runoff_volume = 0.0
    hillslope_volume = 0.0
    outlet_volume = 0.0
    max_step_balance_error = 0.0

    for t in 1:ntime
        depth = view(runoff_depth, :, t)
        map_runoff!(qinst, depth, cell_area, cell_to_reach)
        @assert isapprox(sum(qinst), sum(depth .* cell_area); atol=1e-12)

        route_hillslope!(routed_depth, hillslope, depth)
        map_runoff!(qlat, routed_depth, cell_area, cell_to_reach)
        old_storage = sum(reach_storage(river))
        step!(river, qlat)

        outlet_q[t] = sum(discharge(river)[outlets(net)])
        residual = sum(reach_storage(river)) - old_storage -
                   (sum(qlat) - outlet_q[t]) * dt
        max_step_balance_error = max(max_step_balance_error, abs(residual))

        @assert all(isfinite, discharge(river))
        @assert all(>=(0), discharge(river))
        runoff_volume += sum(qinst) * dt
        hillslope_volume += sum(qlat) * dt
        outlet_volume += outlet_q[t] * dt
    end

    # 5. 验证：坡面汇流守恒；河道出流与末时刻蓄水闭合
    hillslope_error = hillslope_volume - runoff_volume
    river_error = outlet_volume + sum(reach_storage(river)) - hillslope_volume
    tolerance = 1e-8 * runoff_volume
    @assert abs(hillslope_error) <= tolerance
    @assert abs(river_error) <= tolerance
    @assert max_step_balance_error <= tolerance

    @printf("runoff volume:        %.3f m³\n", runoff_volume)
    @printf("peak outlet flow:     %.3f m³/s\n", maximum(outlet_q))
    @printf("hillslope error:      %.3e m³\n", hillslope_error)
    @printf("river balance error:  %.3e m³\n", river_error)
    @printf("max step error:       %.3e m³\n", max_step_balance_error)
    println("validation passed")
end

main()
