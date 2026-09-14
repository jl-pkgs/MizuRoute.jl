function demo_model(method)
    net = RiverNetwork([1, 2, 3], [2, 3, 0])
    p = ReachParameters(3;
        length=[3000.0, 4500.0, 6000.0],
        slope=[0.002, 0.0015, 0.001],
        bottom_width=[8.0, 12.0, 18.0],
        side_slope=1.0,
        floodplain_slope=20.0,
        bankfull_depth=[2.0, 2.5, 3.0],
        mann_n=0.035,
        irf_velocity=1.2,
        irf_diffusivity=300.0,
    )
    RoutingModel(method, net, p; dt=300.0)
end

@testset "Routing methods" begin
    methods = (
        Accumulation(),
        IRF(),
        LagrangianKWT(),
        EulerKinematicWave(cells_per_reach=10),
        MuskingumCunge(),
        DiffusiveWave(cells_per_reach=10),
    )

    for method in methods
        model = demo_model(method)
        for t in 1:100
            qlat = t <= 30 ? [5.0, 0.5, 0.2] : [0.2, 0.1, 0.05]
            vold = copy(reach_storage(model))
            step!(model, qlat)
            @test all(isfinite, discharge(model))
            @test all(>=(0), discharge(model))
            @test all(isfinite, reach_storage(model))
            @test all(>=(0), reach_storage(model))
            for i in model.network.order
                # Because routing is headwater-to-outlet, current upstream
                # discharge is already available when checking reach i.
                qin = sum((discharge(model)[j] for j in upstream_indices(model.network, i)); init=0.0)
                err = water_balance_error(vold[i], reach_storage(model)[i],
                    qin, qlat[i], discharge(model)[i], model.dt)
                @test abs(err) <= 1e-7 * max(1.0, vold[i] + (qin + qlat[i]) * model.dt)
            end
        end
    end
end

@testset "Series API and reset" begin
    model = demo_model(MuskingumCunge())
    qlat = zeros(3, 20)
    qlat[1, :] .= 2.0
    q = route_series(model, qlat)
    @test size(q) == size(qlat)
    @test all(isfinite, q)
    reset!(model)
    @test model.time == 0.0
    @test all(==(0.0), discharge(model))
end
