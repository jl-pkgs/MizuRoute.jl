= 测试、验证与当前边界

== 已编写的自动测试

`test/runtests.jl` 包括以下测试组：

- 河网拓扑：headwater、outlet、直接上游集合、拓扑顺序和环路检测；
- 水力几何：面积--水深互逆、Manning 流量--正常水深互逆；
- 波速和扩散系数的有限性与非负性；
- gamma-UH 权重非负且总和为 1；
- grid/HRU runoff 到 reach lateral inflow 的守恒映射；
- Accumulation、IRF、Lagrangian KWT、Euler KW、MC、DW 六种方案在简单三级河网上的稳定运行；
- 每个时间步的非负流量、非负库容、有限数值；
- 每一河段逐时间步水量守恒；
- `route_series` 与 `reset!` 接口。

运行方式：

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

== 与 mizuRoute 的一致性层级

需要区分三种“一致性”：

+ *方程一致性*：IRF Green function、Muskingum--Cunge 系数以及 Diffusive Wave weighted finite-difference 形式来自当前 mizuRoute 技术说明。该层已按公式实现。
+ *算法结构一致性*：河段按 upstream-to-downstream 顺序计算；侧向入流与上游出流组成当前河段边界；动态水力参数和 sub-stepping 被保留。
+ *逐位一致性*：相同 Cameo testCase、相同参数、相同时间步下，Julia 与 Fortran 的每个 reach / timestep 输出达到指定误差阈值。

第三层目前应视为独立的 regression milestone，而不能仅凭 Julia 自身单元测试宣称完成。尤其 Lagrangian KWT 的 Julia 状态采用 conservative characteristic packets，而原 mizuRoute 继承了 TopNet 的历史 wave-array bookkeeping，因此两者不应预期 bit-for-bit 相同。

== 建议的官方回归流程

当取得 mizuRoute Cameo testCase 的标准输入和 Fortran 输出后，建议为每一种 routing scheme 生成：

$ epsilon_Q(i,t) = Q_("Julia")(i,t) - Q_("Fortran")(i,t), $

并报告：最大绝对误差、相对误差、RMSE、累计体积误差以及 outlet hydrograph 的 NSE/KGE。MC、IRF、Euler KW 和 DW 应优先要求严格数值接近；KWT 可先要求水量守恒、峰现时间和 hydrograph 形状一致，再进一步对齐原 wave bookkeeping。

== Typst 编译

本文档不依赖 Typst Universe 包，可离线编译：

```bash
cd docs
typst compile main.typ MizuRoute.pdf
```

推荐 Typst >= 0.15.1。若需要在论文中使用更多物理量排版宏，可在后续版本中选择性引入 `physica`，但核心技术文档刻意保持零外部依赖。

== GitHub Actions 验证状态

当前仓库已经通过 GitHub Actions 的真实运行测试：Julia 1.10、1.11 和 1.12 三个版本均可成功加载 `MizuRoute`，且 `Pkg.test()` 全部通过；Typst 0.15.1 也可成功编译 `docs/main.typ`，生成的 `MizuRoute.pdf` 由 workflow 作为 artifact 上传。

除运行时 CI 外，仓库还保留了静态源代码检查和独立 Python 数值镜像 sanity check，结果位于 `validation/`。这些内部测试验证了代码自身的一致性，但不能替代与官方 Fortran mizuRoute Cameo/testCase 的逐 reach 数值回归。后者仍是下一阶段最重要的参考验证。
