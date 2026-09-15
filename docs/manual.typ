#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "", header: "")

#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 20pt, weight: "bold")[MizuRoute.jl 使用手册]
  #v(0.4em)
  #text(size: 11pt)[输入数据准备、模型运行与结果读取]
  #v(0.5em)
  #text(size: 9pt)[适用于 MizuRoute.jl v0.1.0]
]

#v(1em)
#outline(title: [目录], depth: 3)
#pagebreak()

= 使用范围

MizuRoute.jl 是供其他 Julia 水文模型调用的轻量河道路由内核。用户负责读取和整理原始数据；本包接收 Julia 数组，不直接解析 CSV、NetCDF、控制文件或重启文件。一次计算包含四步：建立河网、准备河段参数、将产流转换为河段侧向入流、逐时段执行路由。

输入和输出统一采用 SI 单位：时间为 s，长度为 m，面积为 m²，流量为 m³/s。

= 安装与首次运行

== 环境要求

- Julia 1.10、1.11 或 1.12；
- MizuRoute.jl 源码；
- 本包没有第三方 Julia 运行依赖。

在仓库根目录初始化环境：

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

运行自带示例：

```bash
julia --project=. examples/coupling_example.jl
```

运行测试：

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

若从另一个 Julia 项目调用本地源码，可在该项目环境中执行：

```julia
using Pkg
Pkg.develop(path="/path/to/MizuRoute.jl")
```

随后即可 `using MizuRoute`。

= 输入数据准备

所有按河段组织的数组必须使用相同的行顺序；模型输出也沿用该顺序。建议先建立一张“一行一个河段”的主表，再从同一张表生成河网和参数数组。

== 河网拓扑

每个河段需要两个整数：

#table(
  columns: (1.2fr, 2.8fr),
  [字段], [含义],
  `reach_id`, [河段唯一标识，可不连续],
  `downstream_id`, [直接下游河段的 `reach_id`；小于或等于 0 表示流域出口],
)

例如：

```julia
reach_id      = [101, 102, 103, 104]
downstream_id = [103, 103, 104,   0]
net = RiverNetwork(reach_id, downstream_id)
```

该河网表示 101、102 汇入 103，103 汇入 104，104 为出口。构造时会自动生成内部索引和上游至下游的拓扑顺序。

输入必须满足：

- `reach_id` 唯一；
- 每个正的 `downstream_id` 都出现在 `reach_id` 中；
- 河段不能流向自身；
- 河网不能成环；
- 两个数组长度相同且非空。

可用 `headwaters(net)`、`outlets(net)`、`upstream_indices(net, i)` 和 `downstream_index(net, i)` 检查拓扑。它们返回的是 Julia 内部索引，不是 `reach_id`。

== 河段参数

`ReachParameters` 至少需要河长、坡降和河底宽；其他参数可按需要提供。标量会自动扩展到全部河段，向量长度必须等于河段数。

#table(
  columns: (1.6fr, 0.8fr, 2.6fr, 1.3fr),
  [参数], [单位], [含义], [约束/默认值],
  `length`, [m], [河段长度], [> 0，必填],
  `slope`, [-], [河床坡降], [≥ 0，必填],
  `bottom_width`, [m], [主槽底宽], [> 0，必填],
  `side_slope`, [-], [主槽边坡水平/垂直比], [≥ 0；0],
  `floodplain_slope`, [-], [漫滩边坡水平/垂直比], [> 0；1000],
  `bankfull_depth`, [m], [平滩水深], [> 0；1e6],
  `mann_n`, [s·m⁻¹ᐟ³], [Manning 糙率], [> 0；0.035],
  `irf_velocity`, [m/s], [IRF 传播速度], [> 0；1.0],
  `irf_diffusivity`, [m²/s], [IRF 扩散系数], [≥ 0；1000],
)

完整示例：

```julia
p = ReachParameters(length(net);
    length=[4000.0, 3500.0, 8000.0, 12000.0],
    slope=[0.0020, 0.0030, 0.0012, 0.0008],
    bottom_width=[8.0, 6.0, 15.0, 22.0],
    side_slope=1.0,
    floodplain_slope=20.0,
    bankfull_depth=[2.0, 1.8, 2.8, 3.5],
    mann_n=0.035,
    irf_velocity=1.2,
    irf_diffusivity=300.0,
)
```

默认的巨大 `bankfull_depth` 和 `floodplain_slope` 等价于近似关闭漫滩影响。进行实际水力计算时，应尽量提供实测或可靠估算值。

== 侧向入流：直接提供河段流量

若上层模型已输出每个河段或子流域的产流量，整理为长度 `nreach` 的 `qlat`：

```julia
qlat = [2.0, 1.0, 0.5, 0.2]  # m³/s
```

`qlat[i]` 是当前时间步进入第 `i` 个河段的本地侧向入流，不应包含上游河段出流；上游汇流由路由模型自动计算。负侧向入流会被按 0 处理，因此取水应通过 `water_management` 明确传入。

整段序列应组织为 `(nreach, ntime)` 矩阵：每行一个河段，每列一个时间步。

== 侧向入流：由网格或 HRU 径流生成

若上层模型输出径流深率，使用：

```julia
qlat = map_runoff(
    length(net),
    runoff_depth,
    cell_area,
    cell_to_reach;
    weights=fraction,
)
```

各输入含义如下：

#table(
  columns: (1.6fr, 0.8fr, 2.8fr),
  [变量], [单位], [要求],
  `runoff_depth`, [m/s], [每个网格或 HRU 的径流深率],
  `cell_area`, [m²], [对应网格或 HRU 的面积],
  `cell_to_reach`, [-], [接收该单元径流的内部河段索引；0 表示忽略],
  `fraction`, [-], [可选的覆盖或汇入比例；不提供时按 1],
)

计算遵循

$ q_("lat", i) = sum_j R_j A_j W_j. $

注意：`cell_to_reach` 使用 `1:nreach` 的内部索引，而不是任意取值的 `reach_id`。可由 `net.id_to_index[id]` 将河段 ID 转换为内部索引。

常见单位换算：

- mm/时间步：`runoff_depth = runoff_mm / 1000 / dt`；
- mm/h：`runoff_depth = runoff_mm_h / 1000 / 3600`；
- 已为 m³/s：不要再乘面积，直接作为 `qlat`。

反复计算时可预分配以减少内存分配：

```julia
qlat = zeros(length(net))
map_runoff!(qlat, runoff_depth, cell_area, cell_to_reach)
```

== 时间轴与数据检查

路由步长 `dt` 使用秒，且必须与每列输入代表的时段一致。模型不会自动重采样；小时径流通常配 `dt=3600.0`，日径流配 `dt=86400.0`。若产流模块与路由模块步长不同，应在上层模型中先完成守恒聚合或插值。

运行前至少检查：

- 所有输入均为有限数值，不含 `missing`、`NaN` 或 `Inf`；
- 河长、河宽、糙率和面积为正；
- 坡降、径流深率及权重非负；
- 所有河段数组的顺序和长度一致；
- `size(qlat) == (length(net), ntime)`。

= 建模与运行

== 选择路由方法

#table(
  columns: (1.5fr, 3fr),
  [方法], [适用场景],
  `Accumulation()`, [无传播延迟的诊断基线，用于检查拓扑和流量映射],
  `IRF()`, [采用固定速度和扩散系数的脉冲响应路由],
  `LagrangianKWT()`, [Lagrangian 运动波特征线方法],
  `EulerKinematicWave()`, [Euler 隐式运动波方法],
  `MuskingumCunge()`, [一般河道路由的简洁起点，含 CFL 子步进],
  `DiffusiveWave()`, [需要显式表示扩散效应时使用],
)

常用构造方式：

```julia
IRF()
LagrangianKWT(max_packets=20)
EulerKinematicWave(cells_per_reach=20, cfl=1.0)
MuskingumCunge(cfl=0.9)
DiffusiveWave(cells_per_reach=20, alpha=1.0, beta=1.0)
```

Euler 运动波和扩散波要求 `cells_per_reach >= 3`。不确定时可先用 `MuskingumCunge()`，并用 `Accumulation()` 检查输入是否守恒、拓扑是否正确。

== 创建模型

```julia
river = RoutingModel(
    MuskingumCunge(),
    net,
    p;
    dt=3600.0,
    initial_discharge=0.0,
)
```

`initial_discharge` 为所有河段统一的初始流量。默认 `headwater_drain_point=2` 与 mizuRoute 默认语义一致：源头河段的本地径流从河段底部进入，不在该源头河段内经历河道路由。若需从源头河段顶部进入并完整路由，可设为 1。

`active_reaches` 可传入长度为 `nreach` 的布尔向量；非活动上游河段不会计入下游的上游来流。

== 逐时间步运行

逐步调用适合与陆面模型在线耦合：

```julia
Q = Matrix{Float64}(undef, length(net), ntime)

for t in 1:ntime
    qlat = runoff_to_reach(t)
    step!(river, qlat)
    Q[:, t] .= discharge(river)
end
```

`step!` 会修改模型状态，并返回当前各河段出流。`discharge(river)` 返回内部状态向量；需要长期保存时应复制到输出矩阵，不能只保存该向量的引用。

当前状态还可通过下列函数读取：

```julia
discharge(river)     # 河段出流，m³/s
inflow(river)        # 河段上游来流，m³/s
reach_storage(river) # 河段库容，m³
river.time           # 已运行时间，s
```

== 一次运行整段序列

已有 `(nreach, ntime)` 侧向入流矩阵时，可直接执行：

```julia
Q = route_series(river, qlat_series)
```

返回的 `Q` 与输入形状相同，单位为 m³/s。`route_series` 同样会推进模型状态；重复运行同一场景或切换方案前应重置：

```julia
reset!(river)
# 或指定统一初始流量
reset!(river; initial_discharge=1.0)
```

== 坡面延迟

若上层模型只给出瞬时产流，可在河道路由前增加 gamma 单位线：

```julia
hill = GammaUHRouter(
    length(net),
    river.dt;
    shape=2.5,
    scale=7200.0,
)

qlat_delayed = route_hillslope!(hill, qlat_instant)
step!(river, qlat_delayed)
```

`shape` 无量纲，`scale` 单位为 s。若上层模型已表示坡面或地下径流旅行时间，不应再次启用该过程，以免重复延迟。

== 取水与补水

水管理通量使用 mizuRoute 符号约定：正值表示取水，负值表示补水。

```julia
wm = [0.0, 0.2, -0.1, 0.0]  # m³/s
step!(river, qlat; water_management=wm)
```

实际取水受河段可用水量限制。整段计算时，`water_management` 必须与 `qlat` 矩阵同形：

```julia
Q = route_series(river, qlat_series; water_management=wm_series)
```

= 完整可运行示例

将下列内容保存为 `run_routing.jl`，置于仓库根目录后运行。

```julia
using MizuRoute

reach_id = [101, 102, 103, 104]
downstream_id = [103, 103, 104, 0]
net = RiverNetwork(reach_id, downstream_id)

p = ReachParameters(length(net);
    length=[4000.0, 3500.0, 8000.0, 12000.0],
    slope=[0.0020, 0.0030, 0.0012, 0.0008],
    bottom_width=[8.0, 6.0, 15.0, 22.0],
    side_slope=1.0,
    floodplain_slope=20.0,
    bankfull_depth=[2.0, 1.8, 2.8, 3.5],
    mann_n=0.035,
)

river = RoutingModel(MuskingumCunge(), net, p; dt=3600.0)
cell_area = fill(25e6, 8)
cell_to_reach = [1, 1, 2, 2, 3, 3, 4, 4]
Q = Matrix{Float64}(undef, length(net), 48)

for hour in 1:48
    runoff_depth = fill(hour <= 12 ? 1e-6 : 1e-7, 8)
    qlat = map_runoff(length(net), runoff_depth, cell_area, cell_to_reach)
    Q[:, hour] .= step!(river, qlat)
end

println("outlet reach ID = ", net.reach_id[only(outlets(net))])
println("final outlet discharge = ", Q[only(outlets(net)), end], " m³/s")
```

执行：

```bash
julia --project=. run_routing.jl
```

= 结果检查与常见问题

== 建议检查

首次接入新数据时，建议依次检查：

1. 用 `Accumulation()` 运行，确认出口流量量级和河网汇流关系；
2. 检查出口内部索引：`outlets(net)`；
3. 检查 `qlat` 单位是否为 m³/s，尤其避免把 mm/时间步直接当成 m/s；
4. 确认 `cell_to_reach` 使用内部索引而非 `reach_id`；
5. 再切换目标路由方法，检查结果均为有限且非负数；
6. 预热若干时间步，减少零初始状态对分析时段的影响。

单步水量平衡可用 `water_balance_error` 诊断。需要保存旧库容，并对每个河段传入该步的旧库容、新库容、上游来流、侧向入流、出流和 `dt`。数值残差应接近浮点舍入误差。

== 常见报错

#table(
  columns: (2fr, 3fr),
  [现象], [处理],
  [`downstream reach ID ... is not present`], [补齐下游河段，或将出口的下游 ID 设为 0/负数],
  [`river network contains a directed cycle`], [修正河网方向或环路],
  [`parameter length ... != number of reaches`], [使参数向量长度等于 `length(net)`，或改用标量],
  [`cell_area length mismatch`], [使径流、面积、映射和可选权重长度一致],
  [`first dimension must equal number of reaches`], [将时间序列整理为 `(nreach, ntime)`，必要时转置],
  [流量量级异常], [优先检查 mm、m、s、h、day 和面积单位],
)

= 编译本手册

在仓库根目录执行：

```bash
typst compile docs/manual.typ docs/MizuRoute-manual.pdf
```

建议使用 Typst 0.15.1 或更新版本。本手册不依赖 Typst Universe 包，可离线编译。
