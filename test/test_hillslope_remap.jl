@testset "Hillslope UH and runoff remapping" begin
    w = gamma_uh_weights(2.5, 3600.0, 900.0)
    @test all(>=(0), w)
    @test isapprox(sum(w), 1.0; atol=1e-12)
    router = GammaUHRouter(2, 900.0; shape=2.5, scale=3600.0)
    out = route_hillslope!(router, [1.0, 0.0])
    @test all(isfinite, out); @test all(>=(0), out)
    runoff=[1e-6,2e-6,3e-6]; area=[1e6,1e6,2e6]; cell2reach=[1,1,2]
    qlat=map_runoff(2,runoff,area,cell2reach)
    @test isapprox(qlat[1],3.0); @test isapprox(qlat[2],6.0)
end
