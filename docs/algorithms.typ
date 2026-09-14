= 水力学与汇流算法

== 河道断面与 Manning 水力学

主槽采用梯形断面。设水深为 $y$、河底宽度为 $b$、主槽边坡水平/垂直比为 $z_c$，则未超过 bankfull 水深时：

$ A = y (b + z_c y), $
$ B = b + 2 z_c y, $
$ P = b + 2 y sqrt(1 + z_c^2), $
$ R_h = A/P. $

超过 bankfull 后，采用独立的 floodplain 边坡 $z_f$ 扩展面积、顶宽和湿周。均匀流量采用 Manning 公式：

$ Q = 1/n A R_h^(2/3) sqrt(S_0). $

`flow_depth` 使用有界二分法反解给定 $Q$ 的正常水深，因此不依赖外部非线性求解包。

运动波波速按

$ C = partial Q / partial A $

计算，代码中使用围绕当前正常水深的数值微分。扩散波系数采用 mizuRoute 技术说明中的形式

$ D = K^2/(2 Q B), $

结合 $Q=K sqrt(S_0)$ 可写为

$ D = Q/(2 B S_0). $

== Gamma 单位线坡面延迟

mizuRoute 使用 gamma 分布单位线描述产流从坡面到河道入口的延迟 @mizukami2016。连续概率密度为

$ h(t) = 1/(Gamma(k) theta^k) t^(k-1) exp(-t/theta), $

其中 $k$ 为 shape，$theta$ 为 scale。实现中在每个 routing interval 的中点离散采样，并重新归一化，确保

$ sum_m h_m = 1. $

为保持零外部依赖，Gamma 函数使用 Lanczos 近似计算。

== Impulse Response Function (IRF)

IRF 来自一维扩散波方程的 Green 函数 @mizukami2016：

$ Q(x,t) = integral_0^t Q_(in)(t-s) h(x,s) dif s, $

其中

$ h(x,t) = x/(2 t sqrt(pi t D)) exp(-(C t-x)^2/(4 D t)). $

每条河段根据自身长度 $L$、IRF 速度 $C$ 和扩散系数 $D$ 预计算离散 kernel，再对当前及历史输入做卷积。kernel 会重新归一化为单位和，因此不会人为增加或减少水量。

== Lagrangian Kinematic-Wave Tracking (KWT)

原 mizuRoute KWT 源自 TopNet/Goring 的 Lagrangian wave-tracking 思路 @mizukami2016。本 Julia 版本保留其核心物理概念：将一个时间步进入河段的水量表示为沿特征线传播的 conservative packet，其传播速度由当前特征流量的运动波波速决定：

$ dif x / dif t = C(Q). $

每个 packet 保存水量、距出口剩余距离及 characteristic discharge。当剩余距离降至零时，该 packet 在当前时间步成为出口流量。相近 characteristic packets 会自动合并，以限制长模拟中的状态规模。

这一实现是“算法重实现”而不是原 Fortran wave-array 数据结构的逐位复制；因此它满足守恒的 Lagrangian characteristic routing，但不宣称与历史 KWT bookkeeping bit-for-bit 一致。

== Euler Kinematic Wave

忽略扩散项后得到线性化运动波方程：

$ partial Q/partial t + C partial Q/partial x = C q_l. $

当前 mizuRoute 的 Euler-KW 并不是单独的显式迎风求解器，而是令扩散系数 $D=0$，调用与 Diffusive Wave 共用的 `solve_ade` 三对角求解器。Julia 版据此采用相同架构：默认中心差分、全隐式时间权重，上游为给定流量 Dirichlet 边界，下游为保留上一时间层流量梯度的 Neumann 边界。

令

$ C_a = C Delta t / Delta x. $

对于默认全隐式中心差分，内部节点构成三对角系统；代码通过 Thomas 算法求解。虽然全隐式格式不受显式 CFL 稳定性上限约束，Julia 版仍允许依据波速进行 sub-stepping，以便在洪峰快速变化时重新计算水深与运动波波速。

== Muskingum--Cunge

Muskingum--Cunge 采用 mizuRoute 技术说明中的三点显式公式 @cortes2023：

$ O_(t+1) = C_0 I_(t+1) + C_1 I_t + C_2 O_t, $

令 Courant 数为

$ C_n = C Delta t / Delta x, $

Cunge 权重为

$ X = 1/2 (1 - Q/(B S_0 C Delta x)). $

在空间权重 $Y=0.5$ 时：

$ C_0 = (-X + 0.5 C_n)/(1-X+0.5 C_n), $
$ C_1 = ( X + 0.5 C_n)/(1-X+0.5 C_n), $
$ C_2 = (1-X-0.5 C_n)/(1-X+0.5 C_n). $

实现中每一子步根据三点平均流量更新水深、顶宽、波速和 $X$，并通过自动 sub-stepping 控制 Courant 数。

== Diffusive Wave

扩散波控制方程为

$ partial Q/partial t + C partial Q/partial x = D partial^2 Q/partial x^2 + C q_l. $

Julia 版将 mizuRoute `advection_diffusion.f90` 抽成共享的 `advection_diffusion.jl`。令

$ C_a = C Delta t/Delta x, quad C_d = D Delta t/(Delta x)^2. $

`alpha` 对应 mizuRoute 的平流新时间层权重 `wc`，`beta` 对应扩散新时间层权重 `wd`。默认 `alpha=beta=1`，即全隐式。中心差分时，内部节点三对角方程为

$ (-alpha C_a-2 beta C_d) Q_(j-1)^(t+1)
+ (2+4 beta C_d) Q_j^(t+1)
+ (alpha C_a-2 beta C_d) Q_(j+1)^(t+1) = R_j. $

上游采用给定流量 Dirichlet 边界。mizuRoute 当前默认下游边界并非每步强制零梯度，而是保留上一时间层的出口梯度：

$ Q_N^(t+1)-Q_(N-1)^(t+1) = Q_N^t-Q_(N-1)^t. $

三对角系统使用 Thomas 算法求解。当前 mizuRoute 的 `solve_ade` 接口和方程说明包含 lateral flux 参数，但 Fortran RHS 中没有显式加入该项；本 Julia 版遵循其文档化 PDE，默认把 $q_l$ 作为均匀侧向源项加入。该差异在验证章节中单独记录，以便后续做 bug-for-bug 对照。

== 取水与补水

遵循 mizuRoute 约定：正值表示 abstraction，负值表示 injection。正取水依次从：

1. 现有河段储量；
2. 上游入流；
3. 本地侧向入流

中扣除。若请求量超过可用水量，则 `wm_actual` 记录实际可实现的取水量。
