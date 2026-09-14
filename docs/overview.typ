= 设计目标与总体结构

== 目标

MizuRoute.jl 的目标不是逐行复制 mizuRoute 的 Fortran 工程框架，而是将其主要河道汇流算法重新实现为一个可直接嵌入 Julia 水文模型、陆面模式或 PUB 模型的轻量核心库。mizuRoute 的基本工作流是：将水文模型提供的径流深转换并映射为河段侧向入流，必要时先做坡面延迟，再按上游到下游的顺序完成河道路由 @mizukami2016。

本实现保留以下科学计算部分：

- 河网拓扑及上游到下游的拓扑排序；
- 径流深到河段侧向流量的守恒映射；
- gamma 分布单位线坡面汇流；
- IRF、Lagrangian KWT、Euler KW、Muskingum--Cunge 和 Diffusive Wave 五类主要河道路由；
- Manning 水力学、复合梯形断面、波速与扩散系数；
- 河段水量平衡以及可选的取水/补水通量。

以下工程层被有意排除：MPI、PIO、CESM/CTSM coupler、Fortran control-file parser、restart 文件机制及 NetCDF 专用 I/O。这样做可以使 routing kernel 与上层模型彻底解耦。

== 统一耦合接口

上层水文模型只需要提供每个河段在当前时间步的侧向入流 $q_("lat")$，单位为 $"m"^3 "s"^(-1)$。对逐网格模型，可先进行

$ q_("lat", i) = sum_j W_(i,j) R_j A_j, $

其中 $R_j$ 为网格径流深率 $"m" "s"^(-1)$，$A_j$ 为网格面积，$W_(i,j)$ 为网格到河段的守恒映射权重。

运行时接口简化为：

```julia
step!(river, qlat)
Q = discharge(river)
```

因此 routing 模块不需要知道上层产流来自 Richards 方程、VIC、Wflow、BEPS、MarrMot 还是其他模型。

== 河网表示

对每一河段 $i$，仅保存：

- 唯一 reach ID；
- 下游河段索引 $d(i)$；
- 直接上游河段集合 $U(i)$；
- 一个满足所有上游河段先于下游河段的拓扑顺序。

每一时间步对河段 $i$ 的上游流量为

$ Q_("up", i) = sum_(j in U(i)) Q_("out", j). $

由于严格按拓扑顺序遍历，因此在处理 $i$ 时，其所有上游河段的当前时间步出流已经计算完成。

== 水量守恒

所有路由算法最后统一执行河段水量平衡：

$ V^(t+1) = V^t + (Q_("in") + Q_("lat") - Q_("out")) Delta t. $

若候选出流超过当前可用水量，则出流会被限制为可用水量除以 $Delta t$，从而避免负库容。该处理也使不同数值算法拥有一致的质量守恒诊断。
