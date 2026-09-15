= Julia API 与模型对接

== 最小建模流程

```julia
using MizuRoute

net = RiverNetwork(
    [1, 2, 3],
    [2, 3, 0],
)

p = ReachParameters(3;
    length=[3000.0, 4500.0, 6000.0],
    slope=[0.002, 0.0015, 0.001],
    bottom_width=[8.0, 12.0, 18.0],
    side_slope=1.0,
    bankfull_depth=[2.0, 2.5, 3.0],
    mann_n=0.035,
)

river = RoutingModel(
    MuskingumCunge(),
    net,
    p;
    dt=3600.0,
)

for t in eachindex(time)
    qlat = runoff_to_reach(t)
    step!(river, qlat)
    Q[:, t] .= discharge(river)
end
```

== 五种主要方案

```julia
IRF(horizon_factor=6.0)
LagrangianKWT(max_packets=256)
EulerKinematicWave(cells_per_reach=12, cfl=0.85)
MuskingumCunge(cfl=0.9)
DiffusiveWave(cells_per_reach=12, alpha=1.0, beta=1.0)
```

此外 `Accumulation()` 可作为无河道汇流的基线，仅执行上游流量与本地径流的即时累加。

== 参数结构

`ReachParameters` 的核心字段如下：

#table(
  columns: (1.5fr, 1fr, 3fr),
  [参数], [单位], [说明],
  `length`, [m], [河段长度],
  `slope`, [-], [河床坡降],
  `bottom_width`, [m], [主槽底宽],
  `side_slope`, [-], [主槽边坡水平/垂直比],
  `floodplain_slope`, [-], [漫滩边坡水平/垂直比],
  `bankfull_depth`, [m], [bankfull 水深],
  `mann_n`, [-], [Manning 糙率],
  `irf_velocity`, [m/s], [IRF 固定传播速度],
  `irf_diffusivity`, [m²/s], [IRF 固定扩散系数],
)

所有参数既可提供长度为 `nreach` 的向量，也可在构造器中使用标量自动展开。

== 逐网格径流对接

如果上层模型输出网格径流深率 `runoff_depth`，单位为 m/s：

```julia
qlat = map_runoff(
    nreach,
    runoff_depth,
    cell_area,
    cell_to_reach;
    weights=fraction,
)
```

其中 `cell_to_reach[k]` 是第 $k$ 个网格对应的内部河段索引；`weights` 可表示部分覆盖比例。该过程直接完成

$ R ["m/s"] times A ["m"^2] -> Q ["m"^3/"s"]. $

若上层模型已经在 catchment/reach 尺度提供 $"m"^3/"s"$，则无需 remapping，直接将其作为 `qlat` 传给 `step!`。

== 坡面单位线可选耦合

```julia
hill = GammaUHRouter(
    nreach,
    3600.0;
    shape=2.5,
    scale=7200.0,
)

qlat_delayed = route_hillslope!(hill, qlat_instant)
step!(river, qlat_delayed)
```

若上层陆面/水文模型本身已经显式表示了坡面或地下径流旅行时间，应关闭这一层，避免重复延迟。

== 整段时间序列

对形状为 `(nreach, ntime)` 的矩阵：

```julia
Q = route_series(river, qlat)
```

可通过

```julia
reset!(river)
```

清空动态状态后重新运行另一组情景或率定参数。

== 推荐的上层模型边界

建议在上层模型中保持如下职责划分：

$ P -> "ET" -> theta -> R_s, R_b $

由陆面/水文模块计算；然后

$ R_s + R_b -> q_("lat") -> Q_("river") $

完全交给 MizuRoute.jl。这样可以独立替换产流模块和河道模块，也更方便参数率定、EnKF/变分同化及多方案比较。
