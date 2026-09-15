#set document(title: "径流深—坡面汇流—河道汇流案例说明")
#set page(paper: "a4", margin: (x: 2.3cm, y: 2.2cm))
#set text(size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.75em)
#set heading(numbering: "1.")

#let flow-node(title, body) = box(
  width: 84%,
  inset: 8pt,
  radius: 4pt,
  stroke: rgb("4c6a92"),
  fill: rgb("f4f7fb"),
)[
  #strong(title)
  #linebreak()
  #text(size: 9pt, fill: rgb("405066"))[#body]
]

= 径流深—坡面汇流—河道汇流案例说明

本文说明 `runoff_to_river_validation.jl` 的完整计算链：输入网格径流深率，逐网格进行坡面汇流，按河段聚合，沿河网进行扩散波汇流，最后检查水量守恒和数值有效性。

== 总体流程

#figure(
  align(center)[
    #flow-node([网格径流深率], [`depth[k]`，单位 m/s；8 个网格逐时输入])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`route_hillslope!`], [每个网格独立通过 Gamma 单位线，得到延迟后的 `routed_depth[k]`])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`map_runoff!`], [乘网格面积并按 `cell_to_reach` 求和，得到河段侧向入流 `qlat[r]`，单位 m³/s])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([`step!`], [按上游至下游顺序推进河网；本案例采用 Diffusive Wave])
    #v(3pt)
    #text(size: 16pt, fill: rgb("4c6a92"))[↓]
    #v(3pt)
    #flow-node([结果与验证], [出口流量、河道蓄水、坡面与河道水量平衡误差])
  ],
  caption: [案例计算流程],
)

各数组的含义如下：

#table(
  columns: (1.2fr, 1fr, 1fr, 2.5fr),
  inset: 5pt,
  stroke: 0.4pt + rgb("aab4c0"),
  table.header([*变量*], [*长度*], [*单位*], [*含义*]),
  [`depth`], [网格数], [m/s], [当前时刻原始径流深率],
  [`routed_depth`], [网格数], [m/s], [坡面汇流后的网格径流深率],
  [`qinst`], [河段数], [m³/s], [原始径流直接聚合值，仅用于输入水量验证],
  [`qlat`], [河段数], [m³/s], [坡面汇流后进入各河段的侧向流量],
  [`discharge(river)`], [河段数], [m³/s], [当前时刻各河段出口流量],
)

== `route_hillslope!`：逐网格坡面汇流

构造函数为每个网格建立一份独立历史序列，但当前所有网格共用同一组 Gamma 单位线权重：

```julia
hillslope = GammaUHRouter(length(cell_area), dt;
    shape=2.5, scale=3600.0)
```

单位线权重在离散后归一化，使 $sum_j w_j = 1$。每个时刻，`route_hillslope!` 对每个网格执行三步：

+ 删除历史序列最旧值；
+ 将当前径流深率插入序列首端；
+ 计算单位线权重与历史径流的卷积。

对应关系为

$ d_k^(h)(t) = sum_(j=1)^m w_j d_k(t-j+1), $

其中 $d_k(t)$ 是网格 $k$ 的原始径流深率，$d_k^(h)(t)$ 是坡面汇流后的径流深率，$m$ 是单位线长度。案例中的调用为：

```julia
route_hillslope!(routed_depth, hillslope, depth)
```

该函数只描述每个网格产流到达河道前的时间延迟与平滑，不读取 D8/D∞ 流向，也不在网格之间传水，更不负责聚合河段。

== `map_runoff!`：网格到河段的空间聚合

```julia
map_runoff!(qlat, routed_depth, cell_area, cell_to_reach)
```

函数首先将 `qlat` 清零，然后依次处理每个网格。若网格 $k$ 对应河段 $r$，则执行

$ q_("lat",r) += d_k^(h) A_k omega_k, $

其中 $A_k$ 为网格面积，$omega_k$ 为可选权重，默认值为 1。因此单位由 m/s 转为 m³/s。

本案例使用：

```julia
cell_to_reach = [1, 1, 2, 2, 3, 3, 4, 4]
```

即网格 1、2 汇入河段 1，网格 3、4 汇入河段 2，依此类推。`cell_to_reach[k] == 0` 时忽略该网格。

循环中的第一次映射：

```julia
map_runoff!(qinst, depth, cell_area, cell_to_reach)
```

只用于计算坡面汇流前的输入体积和检查空间映射守恒；真正送入河道模型的是坡面汇流后得到的 `qlat`。

== `step!`：河网逐时推进

```julia
step!(river, qlat)
```

`step!` 在一个 `dt` 内完成以下过程：

+ 检查 `qlat` 长度是否等于河段数；
+ 按 `RiverNetwork.order` 从源头向出口遍历；
+ 汇总当前河段所有直接上游河段的出口流量；
+ 读取当前河段的非负侧向入流 `qlat[i]`；
+ 应用可选的取水或补水；
+ 调用所选方法的 `_route_reach!`，更新河段流量剖面、蓄水量和出口流量；
+ 全网完成后将模型时间增加 `dt`。

本案例的河网是

```text
101 ─┐
     ├─> 103 ─> 104 ─> outlet
102 ─┘
```

扩散波方法先依据河段断面、坡度和 Manning 糙率计算水深、波速与扩散系数，再求解一维平流—扩散方程：

$ partial_t Q + c partial_x Q = D partial_x^2 Q. $

空间离散后形成三对角方程组，由 Thomas 算法求解；结果用于更新河段出口流量与蓄水量。当前默认 `headwater_drain_point=2`，因此源头河段 101、102 的本地侧向入流从河段下端进入，不经过该源头河段的扩散波传播；河段 103、104 正常汇流。若希望源头河段也参与河道传播，可构造模型时设置：

```julia
river = RoutingModel(
    DiffusiveWave(cells_per_reach=12), net, p;
    dt=dt, headwater_drain_point=1,
)
```

== 结果验证

案例包含四类检查：

+ *空间映射守恒*：原始网格的 $sum_k d_k A_k$ 等于聚合后的 `sum(qinst)`；
+ *数值有效性*：所有河段流量均有限且非负；
+ *坡面水量守恒*：完整模拟至单位线退水结束后，坡面输出总体积等于输入总体积；
+ *河道水量守恒*：出口累计出流加末时刻河道蓄水等于累计侧向入流。

河网逐时水量残差为

$ epsilon_t = Delta V_t - (sum_r q_("lat",r,t) - sum_o Q_(o,t)) Delta t. $

最终允许误差设为总输入体积的 $10^(-8)$。断言全部通过后输出：

```text
runoff volume:        8640000.000 m³
peak outlet flow:     199.596 m³/s
hillslope error:      1.863e-09 m³
river balance error:  0.000e+00 m³
max step error:       2.083e-10 m³
validation passed
```

== 运行

在项目根目录执行：

```bash
julia --project=. examples/runoff_to_river_validation.jl
typst compile examples/README.typ examples/README.pdf
```
