#set document(title: "Cameo 原始径流—坡面汇流—河道汇流案例")
#set page(paper: "a4", margin: (x: 2.3cm, y: 2.2cm))
#set text(size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.75em)
#set heading(numbering: "1.")

#let flow-node(title, body) = box(
  width: 86%,
  inset: 8pt,
  radius: 4pt,
  stroke: rgb("4c6a92"),
  fill: rgb("f4f7fb"),
)[
  #strong(title)
  #linebreak()
  #text(size: 9pt, fill: rgb("405066"))[#body]
]

= Cameo 原始径流—坡面汇流—河道汇流案例

`runoff_to_river_validation.jl` 直接读取官方 `testCase_cameo_v3.1` 的原始 NetCDF 数据，完成 HRU 径流深输入、逐 HRU 坡面汇流、HRU—河段聚合、真实河网汇流及水量平衡验证。

== 数据来源与规模

案例读取两个原始文件：

+ `input/RUNOFF_case1.nc`：365 天、6745 个 HRU 的瞬时径流深率，单位为 mm/s；
+ `ancillary_data/ntopo_nhdplus_cameo_pfaf.nc`：6736 个有效 HRU 的面积和接收河段，以及 6895 个河段的拓扑、长度和坡度。

径流文件比地形文件多 9 个 HRU。这 9 个 HRU 没有面积和河段映射，不能进入汇流计算；代码按 HRU ID 对齐并保留其余 6736 个有效 HRU。

NetCDF 由全局安装的 `NetCDFTools` 读取：

```julia
using NetCDFTools

runoff_raw = Float64.(
    nc_read(runoff_file, "runoff"; raw=true)
) .* 1e-3 # mm/s → m/s
```

使用 `raw=true` 可直接取得 NetCDF 原始数值，避免整数 ID 的 `_FillValue` 被转换为 `NaN`。

== 总体流程

#figure(
  align(center)[
    #flow-node([原始 Cameo 径流], [`runoff_raw`：6745 HRU × 365 天，mm/s])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([HRU ID 对齐与单位转换], [保留 6736 个有地形映射的 HRU；mm/s 转为 m/s])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`route_hillslope!`], [6736 个 HRU 分别通过 Gamma 单位线，得到 `routed_depth`])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`map_runoff!`], [乘 HRU 面积并按 `hruSegId` 聚合为 6895 个河段的 `qlat`])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`step!`], [按真实河网拓扑逐日推进 IRF 河道汇流])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([结果验证], [出口流量、河道蓄水、逐步及累计水量平衡])
  ],
  caption: [Cameo 原始数据完整汇流流程],
)

== 变量尺寸

本案例有 6745 个原始径流 HRU、6736 个有效地形 HRU、6895 个河段、365 个输入时段和 20 个坡面单位线时段。加入 19 天退水期后共模拟 384 天。

#table(
  columns: (1.35fr, 1.4fr, 0.8fr, 2.5fr),
  inset: 4.5pt,
  stroke: 0.4pt + rgb("aab4c0"),
  table.header([*变量*], [*实际尺寸*], [*单位*], [*含义*]),
  [`runoff_hru`], [`(6745,)`], [—], [径流文件中的 HRU ID],
  [`runoff_raw`], [`(6745, 365)`], [m/s], [单位转换后的全部原始径流深率],
  [`hru_id`], [`(6736,)`], [—], [地形文件中的有效 HRU ID],
  [`cell_area`], [`(6736,)`], [m²], [有效 HRU 面积],
  [`hru_segment`], [`(6736,)`], [—], [每个 HRU 接收河段的真实 ID],
  [`runoff_column`], [`(6736,)`], [—], [有效 HRU 在径流文件中的列索引],
  [`runoff_depth`], [`(6736, 365)`], [m/s], [按 HRU ID 对齐后的径流深率],
  [`reach_id`], [`(6895,)`], [—], [河段真实 ID],
  [`downstream_id`], [`(6895,)`], [—], [下游河段真实 ID],
  [`reach_length`], [`(6895,)`], [m], [河段长度],
  [`reach_slope`], [`(6895,)`], [—], [河段坡度],
  [`cell_to_reach`], [`(6736,)`], [—], [HRU 对应的内部河段索引],
  [`net`], [6895 个河段], [—], [真实河网及拓扑顺序],
  [`hillslope.kernel`], [`(20,)`], [—], [归一化 Gamma 单位线权重],
  [`hillslope.history`], [`(6736,)`，元素 `(20,)`], [m/s], [每个 HRU 独立的径流历史],
  [`depth`], [`(6736,)`], [m/s], [当前时刻原始径流深率],
  [`routed_depth`], [`(6736,)`], [m/s], [当前时刻坡面汇流结果],
  [`qinst`], [`(6895,)`], [m³/s], [原始径流聚合值，仅用于输入验证],
  [`qlat`], [`(6895,)`], [m³/s], [进入各河段的侧向流量],
  [`discharge(river)`], [`(6895,)`], [m³/s], [各河段出口流量],
  [`reach_storage(river)`], [`(6895,)`], [m³], [各河段蓄水量],
  [`outlet_q`], [`(384,)`], [m³/s], [流域出口流量过程],
  [体积及误差变量], [标量], [m³], [累计体积或水量平衡误差],
)

== 原始数据对齐

`NetCDFTools.nc_read` 将径流变量读为“HRU × 时间”矩阵。代码通过真实 HRU ID 建立索引，而不是假定两个文件中的位置完全相同：

```julia
runoff_index = Dict(id => i for (i, id) in pairs(runoff_hru))
runoff_column = [get(runoff_index, id, 0) for id in hru_id]
runoff_depth = runoff_raw[runoff_column, :]
```

随后将 `hruSegId` 转换为 `MizuRoute.jl` 使用的内部河段索引：

```julia
reach_index = Dict(id => i for (i, id) in pairs(reach_id))
cell_to_reach = [get(reach_index, id, 0) for id in hru_segment]
```

== `route_hillslope!`：逐 HRU 坡面汇流

Cameo 参数文件给出 `fshape=2.5`、`tscale=86400 s`。时间步长同为 86400 s：

```julia
hillslope = GammaUHRouter(length(hru_id), dt;
    shape=2.5, scale=86400.0)
```

构造函数生成一组长度为 20 的归一化 Gamma 单位线权重，并为 6736 个 HRU 分别建立独立历史序列。每个时刻，`route_hillslope!`：

+ 删除各 HRU 最旧的历史值；
+ 将当前径流深率插入历史首端；
+ 计算单位线权重与历史径流的卷积。

$ d_k^(h)(t) = sum_(j=1)^20 w_j d_k(t-j+1), quad sum_j w_j=1. $

```julia
route_hillslope!(routed_depth, hillslope, depth)
```

输入和输出均为 m/s。该函数只表示 HRU 产流到达河道前的时间延迟，不在 HRU 之间进行 D8/D∞ 空间传递，也不负责聚合河段。

== `map_runoff!`：HRU 聚合到河段

```julia
map_runoff!(qlat, routed_depth, cell_area, cell_to_reach)
```

函数先将 6895 个河段的 `qlat` 清零，再遍历 6736 个 HRU。若 HRU $k$ 对应河段 $r$，则

$ q_("lat",r) += d_k^(h) A_k. $

径流深率乘面积后，单位由 m/s 转为 m³/s。同一河段接收的多个 HRU 流量直接相加。

循环中另有一次原始径流映射：

```julia
map_runoff!(qinst, depth, cell_area, cell_to_reach)
```

`qinst` 仅用于验证映射前后的瞬时总流量相等；真正进入河道模型的是坡面汇流后的 `qlat`。

== `step!`：真实河网逐日汇流

河网直接由官方 `segId` 和 `downSegId` 构建：

```julia
net = RiverNetwork(reach_id, downstream_id)
river = RoutingModel(
    IRF(legacy_volume_limiter=true), net, p; dt=86400.0,
)
```

`step!(river, qlat)` 在一天内依次执行：

+ 检查 `qlat` 长度是否为 6895；
+ 按 `net.order` 从源头向流域出口遍历；
+ 汇总当前河段所有直接上游河段的出口流量；
+ 加入当前河段侧向流量；
+ 调用 IRF 河段核更新流量历史、河段蓄水和出口流量；
+ 全网完成后将模型时间增加一天。

IRF 以河长、流速 `1.5 m/s` 和扩散系数 `800 m²/s` 构建河段响应函数，对上游来水作卷积延迟。本地侧向流量按 mizuRoute 默认的河段下端入流语义加入。`legacy_volume_limiter=true` 与 Cameo 历史参考程序的限流方式一致。

`ReachParameters` 统一接口要求提供河宽，但 IRF 不使用河宽，因此案例设置 `bottom_width=1.0`，不会参与 IRF 计算。

== 结果验证

案例执行四类检查：

+ *HRU—河段映射守恒*：$sum_k d_k A_k = sum_r q_("inst",r)$；
+ *数值有效性*：所有河段出口流量有限且非负；
+ *坡面水量守恒*：增加 19 天零输入，完整释放坡面单位线历史，使坡面输出总体积等于有效 HRU 输入总体积；
+ *河道水量守恒*：累计出口流量加末时刻河道蓄水等于累计侧向入流。

逐时河网水量残差为

$ epsilon_t = Delta V_t - (sum_r q_("lat",r,t) - sum_o Q_(o,t)) Delta t. $

允许误差设为总输入体积的 $10^(-8)$。当前运行结果为：

```text
raw runoff:           6745 HRUs × 365 days
mapped runoff:        6736 HRUs × 365 days
river network:        6895 reaches
unmapped runoff HRUs: 9
runoff volume:        4073687482.756 m³
peak outlet flow:     902.570 m³/s
hillslope error:      -2.861e-06 m³
river balance error:  2.861e-06 m³
max step error:       3.790e-07 m³
validation passed
```

上述检查验证本实现内部的数值有效性和水量闭合，不表示坡面过程与官方 `dlayRunoff` 逐值一致。官方回归脚本使用参考文件中的 `dlayRunoff`，专门隔离并验证河道汇流核。

== 运行

原始 Cameo 数据默认位于：

```text
cameo/testCase_cameo_v3.1/
```

`NetCDFTools` 使用全局 Julia 环境中的安装，不加入 `MizuRoute.jl` 的核心依赖。运行：

```bash
julia --project=. examples/runoff_to_river_validation.jl
```

也可显式传入数据目录：

```bash
julia --project=. examples/runoff_to_river_validation.jl /path/to/testCase_cameo_v3.1
```

编译本文：

```bash
typst compile examples/README.typ examples/README.pdf
```
