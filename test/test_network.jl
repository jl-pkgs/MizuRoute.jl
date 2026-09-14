@testset "River network" begin
    net = RiverNetwork([10, 20, 30, 40], [30, 30, 40, 0])
    @test length(net) == 4
    @test Set(headwaters(net)) == Set([1, 2])
    @test outlets(net) == [4]
    @test Set(upstream_indices(net, 3)) == Set([1, 2])
    @test downstream_index(net, 3) == 4
    pos = Dict(i => k for (k, i) in pairs(net.order))
    @test pos[1] < pos[3] < pos[4]
    @test pos[2] < pos[3]
    @test_throws ArgumentError RiverNetwork([1, 2], [2, 1])
end
