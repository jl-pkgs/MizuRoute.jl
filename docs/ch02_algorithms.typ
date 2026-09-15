#import "@preview/physica:0.9.8": dd, dv, pdv
#import "@local/modern-cug-report:0.1.3": *

= 水力学与汇流算法

== 河道断面与 Manning 水力学

=== 断面几何

主槽采用梯形断面。设水深为 $y$、河底宽度为 $b$、主槽边坡水平/垂直比为 $z_c$，则未超过 bankfull 水深时：

$ A = y (b + z_c y), $
$ B = b + 2 z_c y, $
$ P = b + 2 y sqrt(1 + z_c^2), $
$ R_h = A/P. $

超过 bankfull 后，采用独立的 floodplain 边坡 $z_f$ 扩展面积、顶宽和湿周。

=== Manning 公式与正常水深

Manning 公式可先写成“河道输水能力”$K$ 与摩阻坡降 $S_f$ 的乘积：

$ K = A/n R_h^(2/3), quad Q = K sqrt(S_f). $

$K$ 将断面面积、湿周和糙率合并为一个量，表示给定断面输送流量的能力；$S_f$ 表示水流克服河床摩阻所需的能量坡降。计算正常水深时采用均匀流假设，即水深沿程不变、水面近似平行于河床，因此 $S_f=S_0$，才得到

$ Q = K sqrt(S_0) = 1/n A R_h^(2/3) sqrt(S_0). $

`flow_depth` 根据该式反解给定 $Q$ 的正常水深。

=== 一维扩散波方程的含义

“一维”是指模型只计算沿河道中心线方向 $x$ 的变化。每个位置只保留断面平均的流量 $Q$、过水面积 $A$ 和水深 $h$，不计算断面内横向和垂向的流速分布。

一维圣维南方程由连续方程和动量方程组成。连续方程表示水量守恒：

$ pdv(A, t) + pdv(Q, x) = q_l. $

动量方程表示水流加速、压力、重力和河床摩阻之间的平衡：

$ pdv(Q, t) + pdv(Q^2/A, x) + g A pdv(h, x) = g A (S_0-S_f) $

$ 1/A pdv(Q, t) + 1/A pdv(Q^2/A, x) + g pdv(h, x) - g S_0 + g S_f = 0 $

从左到右依次是局地惯性、对流惯性、水面压力、重力和摩阻。三种波动近似都保留连续方程，区别只在于动量方程保留哪些项：

#figure(
  caption: [
    三种波动近似在动量方程中的保留项。
  ],
  table(
    columns: (1.2fr, 1.4fr, 1.4fr, 1.4fr, 1.0fr, 1.0fr, 3.0fr),
    rows: (0.9cm, 0.9cm),
    align: center + horizon,
    [近似], [$1/A pdv(Q, t)$], [$1/A pdv(Q^2/A, x)$], [$g pdv(h, x)$], [$-g S_0$], [$+g S_f$], [化简结果],
    [物理项], [局地惯性], [对流惯性], [水面压力], [重力], [摩阻], [],
    [动力波], [✓], [✓], [✓], [✓], [✓], [不作简化，求解完整圣维南方程],
    [扩散波], [×], [×], [✓], [✓], [✓], [$S_f=S_0-pdv(h, x)$],
    [运动波], [×], [×], [×], [✓], [✓], [$S_f=S_0$，$Q=K sqrt(S_0)$],
  ),
) <table_>
#table-note[
   ✓ 表示保留；× 表示忽略，“忽略”是指该项相对于其余项足够小。
]

“扩散波”是指忽略两个惯性项、但保留水面压力梯度的圣维南近似。它描述的不是水分子的扩散，而是高水位河段对相邻低水位河段的推动，使洪水过程在传播中削峰、展宽。

河道中的洪水波会同时发生三件事：以一定速度向下游移动、在移动过程中逐渐展宽，以及接收沿岸坡面来水。将一维扩散波方程在当前水力状态附近线性化，可写成：

$ pdv(Q, t) + C pdv(Q, x) = D pdv(Q, x, 2) + C q_l. $

可将它直接读作

```text
流量随时间变化 + 洪水波向下游移动 = 波形展宽削峰 + 坡面入流
```

各项含义如下：

#table(
  columns: (1.5fr, 2.2fr, 2.4fr),
  [项], [数学含义], [水文含义],
  $pdv(Q, t)$, [固定位置处流量的时间变化率], [描述该位置正在涨水还是退水],
  $C pdv(Q, x)$, [流量沿河道的空间梯度乘波速], [把洪水过程以速度 $C$ 搬向下游],
  $D pdv(Q, x, 2)$, [流量沿程曲线的弯曲程度乘扩散系数], [压低尖锐洪峰并填平低谷],
  $C q_l$, [单位河长侧向入流造成的流量变化], [表示坡面、支沟等沿河道补水],
)

其中 $x$ 是沿河道向下游的距离，$t$ 是时间，$Q(x,t)$ 是流量，$q_l$ 是单位河长的侧向入流。$C$ 的单位为 m/s，决定洪峰何时到达；$D$ 的单位为 m²/s，决定洪峰展宽多快。

最简单的理解是：若 $D=0$ 且 $q_l=0$，一个洪水过程只以速度 $C$ 向下游平移，形状不变；若 $D>0$，尖峰处的二阶导数为负，扩散项会降低峰值，而低谷处的二阶导数为正，扩散项会抬高低谷。因此 $D$ 产生“削峰、展宽”的效果。

这里“线性化”是指：真实的 $C$ 和 $D$ 会随流量和水深变化，但在一个求解步内先根据代表性流量算出，再把它们当作已知系数求解。

=== 波速与扩散系数的推导

下面从 Manning 关系逐步推导 $C$ 和 $D$。这里的 $dd(Q)$ 不是时间导数，而是：当水流状态发生一个很小的变化时，流量随之产生的微小变化。

在河道断面形状和糙率固定时，输水能力 $K$ 是过水面积 $A$ 的函数：

$ K = K(A). $

Manning 关系因此应完整写成一个二元函数：

$ Q(A, S_f) = K(A) sqrt(S_f). $

也就是说，$Q$ 会因 $A$ 改变而改变，也会因 $S_f$ 改变而改变。二元函数的全微分公式是

$ dd(Q) = (pdv(Q, A))_(S_f) dd(A) + (pdv(Q, S_f))_A dd(S_f). $

第一项表示“保持摩阻坡降 $S_f$ 不变，只改变过水面积”；第二项表示“保持面积 $A$ 不变，只改变摩阻坡降”。下面分别计算这两个偏导数。

对 $A$ 求偏导时，$sqrt(S_f)$ 是常数：

$ (pdv(Q, A))_(S_f) = sqrt(S_f) dv(K, A) equiv C. $

这一定义给出波速 $C$。对 $S_f$ 求偏导时，$K(A)$ 是常数：

$ (pdv(Q, S_f))_A = K(A) pdv(sqrt(S_f), S_f) = K/(2 sqrt(S_f)). $

再由 $Q=K sqrt(S_f)$ 得到 $sqrt(S_f)=Q/K$，所以

$ K/(2 sqrt(S_f)) = K/(2 Q/K) = K^2/(2 Q). $

于是流量全微分完整写为

$ dd(Q) = C dd(A) + K^2/(2 Q) dd(S_f). $

接下来把摩阻坡降的变化转换为面积的空间变化。忽略完整动量方程中的惯性项后，

$ S_f = S_0 - pdv(h, x). $

河床坡降 $S_0$ 固定，因此两边取微小变化：

$ dd(S_f) = -dd(pdv(h, x)) = -pdv(dd(h), x). $

水面升高一个微小量 $dd(h)$ 时，增加的过水面积近似等于水面宽度乘水深增量：

$dd(A) = B dd(h)$，故$dd(h) = dd(A)/B$。

在当前水力状态附近进行局部线性化，把 $B$ 暂时视为常数，则

$ dd(S_f) approx -1/B pdv(dd(A), x). $

代回流量全微分：

$ dd(Q) = C dd(A) - K^2/(2 Q B) pdv(dd(A), x). $

定义面积空间梯度前的系数为水动力扩散系数：

$ D equiv K^2/(2 Q B). $

于是

$ dd(Q) = C dd(A) - D pdv(dd(A), x). $

最后说明它如何得到前述平流--扩散方程。将上述微小变化理解为一个微小时间段内的变化，并除以该时间增量：

$ pdv(Q, t) = C pdv(A, t) - D pdv(pdv(A, t), x). $

连续方程给出

$ pdv(A, t) = q_l - pdv(Q, x). $

代入上式：

$ pdv(Q, t) = C (q_l - pdv(Q, x)) - D (pdv(q_l, x) - pdv(Q, x, 2)). $

若一个河段内的侧向入流近似均匀，即 $pdv(q_l,x)=0$，移项后得到

$ pdv(Q, t) + C pdv(Q, x) = D pdv(Q, x, 2) + C q_l. $

$D$ 的量纲为 m²/s。直观上，某处水位高于相邻河段时，水面坡度产生额外压力梯度，使洪水波由高水位区向两侧展宽；$D$ 就衡量这种展宽速度。

由 $S_f=(Q/K)^2$ 还可得

$ D = Q/(2 B S_f). $

在正常水深的均匀流参考状态下 $S_f=S_0$，因此进一步写成

$ D = Q/(2 B S_0). $

该推导使用局部棱柱形河道和缓变水流近似，并在当前水力状态附近线性化。$D$ 并非所有方法都要求用户输入：`DiffusiveWave` 根据当前流量、断面和坡降动态计算它；`EulerKinematicWave` 令 $D=0$；`IRF` 则使用用户提供的固定 `irf_diffusivity`。两者单位和物理含义相同，但计算方式不同。


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

$ dv(x, t) = C(Q). $

每个 packet 保存水量、距出口剩余距离及 characteristic discharge。当剩余距离降至零时，该 packet 在当前时间步成为出口流量。相近 characteristic packets 会自动合并，以限制长模拟中的状态规模。

这一实现是“算法重实现”而不是原 Fortran wave-array 数据结构的逐位复制；因此它满足守恒的 Lagrangian characteristic routing，但不宣称与历史 KWT bookkeeping bit-for-bit 一致。

== Euler Kinematic Wave

忽略扩散项后得到线性化运动波方程：

$ pdv(Q, t) + C pdv(Q, x) = C q_l. $

当前 mizuRoute 的 Euler-KW 并不是单独的显式迎风求解器，而是令扩散系数 $D=0$，调用与 Diffusive Wave 共用的 `solve_ade` 三对角求解器。Julia 版据此采用相同架构：默认中心差分、全隐式时间权重，上游为给定流量 Dirichlet 边界，下游为保留上一时间层流量梯度的 Neumann 边界。

令

$ C_a = C Delta t / Delta x. $

对于默认全隐式中心差分，内部节点构成三对角系统；代码通过 Thomas 算法求解。虽然全隐式格式不受显式 CFL 稳定性上限约束，Julia 版仍允许依据波速进行 sub-stepping，以便在洪峰快速变化时重新计算水深与运动波波速。

== Muskingum--Cunge

Muskingum--Cunge 采用 mizuRoute 技术说明中的三点显式公式 @cortes2023：

$ O_(t+1) = C_0 I_(t+1) + C_1 I_t + C_2 O_t, $

令 Courant 数为

$ C_n = C (Delta t) / (Delta x), $

Cunge 权重为

$ X = 1/2 (1 - Q/(B S_0 C Delta x)). $

在空间权重 $Y=0.5$ 时：

$ C_0 = (-X + 0.5 C_n)/(1-X+0.5 C_n), $
$ C_1 = ( X + 0.5 C_n)/(1-X+0.5 C_n), $
$ C_2 = (1-X-0.5 C_n)/(1-X+0.5 C_n). $

实现中每一子步根据三点平均流量更新水深、顶宽、波速和 $X$，并通过自动 sub-stepping 控制 Courant 数。

== Diffusive Wave

扩散波控制方程为

$ pdv(Q, t) + C pdv(Q, x) = D pdv(Q, x, 2) + C q_l. $

Julia 版将 mizuRoute `advection_diffusion.f90` 抽成共享的 `advection_diffusion.jl`。令

$ C_a = C (Delta t) / (Delta x), quad C_d = D (Delta t) / (Delta x)^2. $

`alpha` 对应 mizuRoute 的平流新时间层权重 `wc`，`beta` 对应扩散新时间层权重 `wd`。默认 `alpha=beta=1`，即全隐式。中心差分时，内部节点三对角方程为

$ (-alpha C_a-2 beta C_d) Q_(j-1)^(t+1)
+ (2+4 beta C_d) Q_j^(t+1)
+ (alpha C_a-2 beta C_d) Q_(j+1)^(t+1) = R_j. $

上游采用给定流量 Dirichlet 边界。mizuRoute 当前默认下游边界并非每步强制零梯度，而是保留上一时间层的出口梯度：

$ Q_N^(t+1)-Q_(N-1)^(t+1) = Q_N^t-Q_(N-1)^t. $

三对角系统使用 Thomas 算法求解。需要注意，上式中的 $C q_l$ 是理论方程中的沿程源项；当前 mizuRoute 的 `solve_ade` 并未将它写入三对角方程右端，而是在河道路由完成后把该河段的总侧向入流加到出口流量。本 Julia 版为保持回归一致也采用这一处理。

== 取水与补水

遵循 mizuRoute 约定：正值表示 abstraction，负值表示 injection。正取水依次从：

1. 现有河段储量；
2. 上游入流；
3. 本地侧向入流

中扣除。若请求量超过可用水量，则 `wm_actual` 记录实际可实现的取水量。
