@testset "Shared ADE solver" begin
    qconst = fill(3.5, 8)
    q1 = MizuRoute._solve_ade(qconst, 5000.0, 300.0, 3.5, 1.2, 0.0, 0.0)
    @test all(isapprox.(q1, 3.5; rtol=1e-12, atol=1e-12))

    qzero = zeros(8)
    q2 = MizuRoute._solve_ade(qzero, 5000.0, 300.0, 0.0, 1.2, 200.0, 0.0)
    @test all(isapprox.(q2, 0.0; atol=1e-14))

    q3 = MizuRoute._solve_ade(qzero, 5000.0, 300.0, 2.0, 1.2, 200.0, 0.0)
    @test all(isfinite, q3)
    @test all(>=(0.0), q3)
    @test q3[1] ≈ 2.0

    q4 = MizuRoute._solve_ade(qzero, 5000.0, 300.0, 0.0, 1.2, 200.0, 1.0)
    @test maximum(q4) > 0.0
end
