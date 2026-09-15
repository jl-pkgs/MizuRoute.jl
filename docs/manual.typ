#import "@local/modern-cug-report:0.1.3": *
#import "@preview/physica:0.9.8": dd, pdv
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

MizuRoute.jl 是供其他 Julia 水文模型调用的轻量河道汇流内核。用户负责读取和整理原始数据；本包接收 Julia 数组，不直接解析 CSV、NetCDF、控制文件或重启文件。一次计算包含四步：建立河网、准备河段参数、将产流转换为河段侧向入流、逐时段执行汇流。

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
  `irf_velocity`, [m/s], [IRF 河道响应传播速度], [> 0；1.0],
  `irf_diffusivity`, [m²/s], [IRF 河道响应纵向扩散系数], [≥ 0；1000],
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

=== `irf_velocity` 与 `irf_diffusivity`

这两个参数仅供 `IRF()` 使用。IRF 将河段中的流量传播近似为一维平流--扩散过程：

$ pdv(Q, t) + C pdv(Q, x) = D pdv(Q, x, 2), $

其中 $C$ 即 `irf_velocity`，$D$ 即 `irf_diffusivity`。对长度为 $L$ 的河段，其单位脉冲响应为

$ h(L, t) = L/(2 t sqrt(pi D t)) exp(-(C t - L)^2/(4 D t)), quad t > 0. $

上游入流通过卷积变为河段出流：

$ Q_("out")(t) = integral_0^t Q_("in")(t-s) h(L, s) dd(s). $

因此两个参数的作用可直接由响应时间近似量化：

$ t_("mean") = L/C, quad sigma_t = sqrt(2 D L/C^3). $

- `irf_velocity`（$C$）是等效洪水波传播速度，并非断面平均水流速度。$C$ 越大，平均到达时间 $L/C$ 越短，洪峰越早出现。
- `irf_diffusivity`（$D$）是河段尺度的等效纵向扩散系数。$D$ 越大，响应时间标准差 $sigma_t$ 越大，过程线越宽、洪峰越低；$D$ 越小，响应越集中。

例如 $L=8000$ m、$C=1.2$ m/s、$D=300$ m²/s 时，平均传播时间约为 1.85 h，响应时间标准差约为 0.46 h。率定时可先用洪峰到达时间约束 $C$，再用过程线宽度和峰值约束 $D$。`IRF()` 以外的汇流方法不读取这两个参数。

默认的巨大 `bankfull_depth` 和 `floodplain_slope` 等价于近似关闭漫滩影响。进行实际水力计算时，应尽量提供实测或可靠估算值。

== 侧向入流：直接提供河段流量

若上层模型已输出每个河段或子流域的产流量，整理为长度 `nreach` 的 `qlat`：

```julia
qlat = [2.0, 1.0, 0.5, 0.2]  # m³/s
```

`qlat[i]` 是当前时间步进入第 `i` 个河段的本地侧向入流，不应包含上游河段出流；上游汇流由汇流模型自动计算。负侧向入流会被按 0 处理，因此取水应通过 `water_management` 明确传入。

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

汇流步长 `dt` 使用秒，且必须与每列输入代表的时段一致。模型不会自动重采样；小时径流通常配 `dt=3600.0`，日径流配 `dt=86400.0`。若产流模块与汇流模块步长不同，应在上层模型中先完成守恒聚合或插值。

运行前至少检查：

- 所有输入均为有限数值，不含 `missing`、`NaN` 或 `Inf`；
- 河长、河宽、糙率和面积为正；
- 坡降、径流深率及权重非负；
- 所有河段数组的顺序和长度一致；
- `size(qlat) == (length(net), ntime)`。

= 当前坡面汇流处理 <hillslope-routing>

== 作用与计算顺序

分布式水文模型通常先计算各网格产流，再计算水从坡面进入河道的时间过程。MizuRoute.jl 将这两个步骤明确分开：`map_runoff!` 只做空间汇总，`route_hillslope!` 再做时间延迟，最后由 `step!` 完成河道路由。

```text
网格产流深率
  → 按面积汇总为河段瞬时入流
  → Gamma 单位线延迟
  → 河段侧向入流
  → 河道路由
```

对第 $n$ 个时间步，网格径流首先转换为第 $i$ 个河段的瞬时入流：

$ q_("inst", i)^n = sum_(j: r_j=i) R_j^n A_j W_j, $

其中 $R_j$ 为网格径流深率，$A_j$ 为网格面积，$W_j$ 为可选权重，$r_j$ 为 `cell_to_reach[j]`。这一步只有单位转换和求和，不包含坡面传播时间。

== Gamma 单位线

坡面传播采用 Gamma 概率密度函数：

$ f(t; a, theta) = 1/(Gamma(a) theta^a) t^(a-1) exp(-t/theta), quad t > 0, $

其中 `shape` 对应形状参数 $a$，`scale` 对应时间尺度 $theta$，单位为 s。其时间特征为

$ t_("mean") = a theta, quad t_("peak") = (a-1) theta, quad sigma_t = sqrt(a) theta, $

其中峰值时间公式适用于 $a>1$。因此，增大 `scale` 会整体推迟并展宽入河过程；增大 `shape` 也会延长平均传播时间，并改变过程线形状。

当前实现不是直接积分 Gamma 分布，而是在每个路由时段的中点采样：

$ t_k = (k-1/2) Delta t, quad w_k = (f(t_k) Delta t)/(sum_(m=1)^N f(t_m) Delta t). $

默认取

$ N = max(8, ceil(8 a theta / Delta t)). $

权重会重新归一化，满足 $sum_k w_k=1$。因此，只要继续运行到单位线尾部全部流出，坡面延迟不会改变总水量。

每个河段分别保存自身的历史入流，但当前所有河段共用同一组 Gamma 权重。第 $n$ 步进入河道的侧向流量为

$ q_("lat", i)^n = sum_(k=1)^N w_k q_("inst", i)^(n-k+1). $

== 运行方法

建议预分配瞬时入流和延迟入流：

```julia
hill = GammaUHRouter(
    length(net),
    river.dt;
    shape=2.5,
    scale=7200.0,
)

qinst = zeros(length(net))
qlat = zeros(length(net))

for t in 1:ntime
    map_runoff!(qinst, runoff[:, t], cell_area, cell_to_reach)
    route_hillslope!(qlat, hill, qinst)
    Q[:, t] .= step!(river, qlat)
end
```

`GammaUHRouter` 是有记忆的动态状态。重新运行情景时，除执行 `reset!(river)` 外，还需重新构造 `hill`。模拟结束时若单位线尾部仍有水量保存在历史状态中，应继续输入零径流运行若干步，或将该部分作为末时刻未排出的坡面储量处理。

== flowdir 在哪里使用

当前模型不接收 flowdir，也不在运行时沿 D8/D∞ 方向逐格追踪水流。flowdir 应在预处理中用于确定 `cell_to_reach`，即每个网格最终汇入哪个河段。运行时只读取这一映射结果。

当前顺序是“先按河段汇总网格径流，再使用统一 Gamma 单位线”。因此，网格到河道的具体距离会在汇总时丢失：近河网格和远河网格使用相同的延迟分布。Gamma 单位线只能统计性地表示一个河段汇水区的综合入河时间，不能显式区分每个网格的旅行时间。

原 Fortran mizuRoute 采用相同的总体思路：由 `hruSegId` 指定 HRU 汇入的河段，再用空间统一的 Gamma 单位线延迟；它同样不在运行期执行 flowdir 栅格汇流。不同之处是原版通过 Gamma 累积分布函数计算每个时间段的概率，而当前 Julia 实现使用中点采样。

若要显式考虑网格距离，应在求和前根据 `distance_to_channel` 或 `travel_time` 分别延迟各网格径流，再汇总到河段。该改进记录于仓库根目录的 `TODO.md`。

== 何时跳过坡面汇流

若上层水文模型输出的是未经传播的网格产流，应执行 `map_runoff!` 和 `route_hillslope!`。若上层模型已经沿坡面或地下路径计算了入河时间，则应直接把入河流量作为 `qlat` 传给 `step!`，避免重复延迟。

= 建模与运行

== 选择汇流方法

不同方法使用不同的河段参数。`RoutingModel` 统一接收一个完整的 `ReachParameters`，但所选方法只读取下表中的字段；未使用的字段可保留默认值。

#table(
  columns: (1.4fr, 2.1fr, 2.7fr),
  [方法], [使用的河段参数], [适用场景],
  `Accumulation()`, [无], [无传播延迟的诊断基线，用于检查拓扑和流量映射],
  `IRF()`, [`length`、`irf_velocity`、`irf_diffusivity`], [采用固定速度和扩散系数的脉冲响应汇流],
  `LagrangianKWT()`, [`length`、`slope`、`mann_n`], [Lagrangian 运动波特征线方法],
  `EulerKinematicWave()`, [`length`、`slope`、断面参数、`mann_n`], [Euler 隐式运动波方法],
  `MuskingumCunge()`, [`length`、`slope`、断面参数、`mann_n`], [一般河道汇流的简洁起点，含 CFL 子步进],
  `DiffusiveWave()`, [`length`、`slope`、断面参数、`mann_n`], [需要显式表示扩散效应时使用],
)

表中“断面参数”包括 `bottom_width`、`side_slope`、`floodplain_slope` 和 `bankfull_depth`。方法构造器中的 `cfl`、`cells_per_reach`、`alpha`、`beta`、`max_packets` 等是数值算法设置，不属于河段物理参数。

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

`initial_discharge` 为所有河段统一的初始流量。默认 `headwater_drain_point=2` 与 mizuRoute 默认语义一致：源头河段的本地径流从河段底部进入，不在该源头河段内经历河道汇流。若需从源头河段顶部进入并完整汇流，可设为 1。

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
5. 再切换目标汇流方法，检查结果均为有限且非负数；
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
