using MizuRoute

net = RiverNetwork([101,102,103,104],[103,103,104,0])
p = ReachParameters(length(net);
    length=[4000.0,3500.0,8000.0,12000.0],
    slope=[0.0020,0.0030,0.0012,0.0008],
    bottom_width=[8.0,6.0,15.0,22.0],
    side_slope=1.0,
    floodplain_slope=20.0,
    bankfull_depth=[2.0,1.8,2.8,3.5],
    mann_n=0.035,
    irf_velocity=1.2,
    irf_diffusivity=300.0,
)
river = RoutingModel(DiffusiveWave(cells_per_reach=12),net,p;dt=3600.0)
cell_area=fill(25e6,8)
cell_to_reach=[1,1,2,2,3,3,4,4]
for hour in 1:48
    runoff_depth=hour<=12 ? fill(1e-6,8) : fill(1e-7,8)
    qlat=map_runoff(length(net),runoff_depth,cell_area,cell_to_reach)
    step!(river,qlat)
    @show hour discharge(river)
end
