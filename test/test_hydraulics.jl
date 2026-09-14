@testset "Hydraulics" begin
    b,z,zf,bf=8.0,1.5,20.0,2.0
    for y in (0.1,1.0,2.0,3.5)
        a=flow_area(y,b,z;zf=zf,bankfull_depth=bf)
        y2=water_height(a,b,z;zf=zf,bankfull_depth=bf)
        @test isapprox(y,y2;rtol=1e-10,atol=1e-10)
        @test top_width(y,b,z;zf=zf,bankfull_depth=bf)>0
        @test wetted_perimeter(y,b,z;zf=zf,bankfull_depth=bf)>0
    end
    q=35.0; slope=1e-3; n=0.035
    y=flow_depth(q,b,z,slope,n;zf=zf,bankfull_depth=bf)
    q2=manning_discharge(y,b,z,slope,n;zf=zf,bankfull_depth=bf)
    @test isapprox(q,q2;rtol=1e-6)
    @test celerity(q,b,z,slope,n;zf=zf,bankfull_depth=bf)>0
    @test diffusivity(q,b,z,slope,n;zf=zf,bankfull_depth=bf)>=0
end
