using MizuRoute
using NetCDFTools
using Printf

function main(root=isempty(ARGS) ?
    normpath(joinpath(@__DIR__, "..", "cameo", "testCase_cameo_v3.1")) : ARGS[1])
    runoff_file = joinpath(root, "input", "RUNOFF_case1.nc")
    topo_file = joinpath(root, "ancillary_data", "ntopo_nhdplus_cameo_pfaf.nc")
    isfile(runoff_file) || error("Cameo runoff file not found: $runoff_file")
    isfile(topo_file) || error("Cameo topology file not found: $topo_file")

    # 1. 读取 Cameo 原始 HRU 径流深率和河网数据
    runoff_hru = Int.(nc_read(runoff_file, "hru"; raw=true))
    runoff_raw = Float64.(nc_read(runoff_file, "runoff"; raw=true)) .* 1e-3 # mm/s → m/s
    hru_id = Int.(nc_read(topo_file, "HRUid"; raw=true))
    cell_area = Float64.(nc_read(topo_file, "area"; raw=true))
    hru_segment = Int.(nc_read(topo_file, "hruSegId"; raw=true))
    reach_id = Int.(nc_read(topo_file, "segId"; raw=true))
    downstream_id = Int.(nc_read(topo_file, "downSegId"; raw=true))
    reach_length = Float64.(nc_read(topo_file, "length"; raw=true))
    reach_slope = Float64.(nc_read(topo_file, "slope"; raw=true))

    # 径流文件多出 9 个无地形映射的 HRU；按 HRU ID 对齐有效输入
    runoff_index = Dict(id => i for (i, id) in pairs(runoff_hru))
    runoff_column = [get(runoff_index, id, 0) for id in hru_id]
    all(>(0), runoff_column) || error("topology contains HRUs absent from runoff input")
    runoff_depth = runoff_raw[runoff_column, :]

    reach_index = Dict(id => i for (i, id) in pairs(reach_id))
    cell_to_reach = [get(reach_index, id, 0) for id in hru_segment]
    all(>(0), cell_to_reach) || error("an HRU points to an unknown reach")

    # 2. 构建 6895 段真实河网；IRF 使用 Cameo 参数
    net = RiverNetwork(reach_id, downstream_id)
    p = ReachParameters(length(net);
        length=reach_length,
        slope=max.(reach_slope, 1e-6),
        bottom_width=1.0, # IRF 不使用河宽，仅满足统一参数接口
        irf_velocity=1.5,
        irf_diffusivity=800.0,
    )
    dt = 86400.0

    # 3. 6736 个有效 HRU 分别进行 Gamma 单位线坡面汇流
    hillslope = GammaUHRouter(length(hru_id), dt; shape=2.5, scale=86400.0)
    input_steps = size(runoff_depth, 2)
    ntime = input_steps + length(hillslope.kernel) - 1

    # 4. 坡面出流聚合到河段，再进行 IRF 河道汇流
    river = RoutingModel(IRF(legacy_volume_limiter=true), net, p; dt=dt)
    zero_depth = zeros(length(hru_id))
    routed_depth = similar(zero_depth)
    qinst = zeros(length(net))
    qlat = similar(qinst)
    outlet_q = zeros(ntime)
    runoff_volume = 0.0
    hillslope_volume = 0.0
    outlet_volume = 0.0
    max_step_balance_error = 0.0

    for t in 1:ntime
        depth = t <= input_steps ? view(runoff_depth, :, t) : zero_depth
        map_runoff!(qinst, depth, cell_area, cell_to_reach)
        @assert isapprox(sum(qinst), sum(depth .* cell_area); rtol=1e-12, atol=1e-10)

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

    # 5. 验证坡面和河道水量守恒
    hillslope_error = hillslope_volume - runoff_volume
    river_error = outlet_volume + sum(reach_storage(river)) - hillslope_volume
    tolerance = 1e-8 * runoff_volume
    @assert abs(hillslope_error) <= tolerance
    @assert abs(river_error) <= tolerance
    @assert max_step_balance_error <= tolerance

    @printf("raw runoff:           %d HRUs × %d days\n", size(runoff_raw)...)
    @printf("mapped runoff:        %d HRUs × %d days\n", size(runoff_depth)...)
    @printf("river network:        %d reaches\n", length(net))
    @printf("unmapped runoff HRUs: %d\n", length(runoff_hru) - length(hru_id))
    @printf("runoff volume:        %.3f m³\n", runoff_volume)
    @printf("peak outlet flow:     %.3f m³/s\n", maximum(outlet_q))
    @printf("hillslope error:      %.3e m³\n", hillslope_error)
    @printf("river balance error:  %.3e m³\n", river_error)
    @printf("max step error:       %.3e m³\n", max_step_balance_error)
    println("validation passed")
end

main()
