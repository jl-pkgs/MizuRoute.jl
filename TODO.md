# TODO

## 考虑网格到河道距离的坡面汇流

### 当前问题

当前流程为：

```text
网格径流 → map_runoff 按河段求和 → GammaUHRouter 统一延迟 → 河道汇流
```

网格求和后，其到河道的距离信息已经丢失。现有 `GammaUHRouter` 对所有河段使用同一组 `shape` 和 `scale`，只能笼统表示坡面汇流延迟，不能区分近河网格与远河网格。

### 原 Fortran 的处理

原 mizuRoute 同样不在运行时沿 D8/D∞ flowdir 逐格汇流。flowdir 仅在预处理中用于生成：

- `hruSegId`：HRU 汇入的河段；
- `downSegId`：河段下游关系；
- 网格或水文模型 HRU 到河网 HRU 的面积权重。

运行时先将径流映射到河网 HRU，再用空间统一的 Gamma 单位线延迟，最后执行河道汇流。因此原版也没有显式考虑每个网格到河道距离的差异。

相关源码：

- `../MizuRoute/route/build/src/main_route.f90`
- `../MizuRoute/route/build/src/process_remap.f90`
- `../MizuRoute/route/build/src/process_param.f90`
- `../MizuRoute/route/build/src/basinUH.f90`

### 改进方案

flowdir 仍只用于预处理，不放入时间循环。预先为每个网格计算并保存：

- 汇入河段索引 `cell_to_reach`；
- 到河道的流路长度 `distance_to_channel`，或旅行时间 `travel_time`；
- 网格有效面积或重映射权重。

运行时应先按网格旅行时间延迟，再汇总为河段侧向入流：

```text
网格径流
  → 按 distance_to_channel / travel_time 延迟
  → 汇总为各河段 qlat
  → 河道汇流
```

可令网格平均旅行时间为

```text
τ[k] = distance_to_channel[k] / hillslope_velocity[k]
```

并以 `τ[k]` 构造网格单位线；若采用 Gamma 分布，可固定 `shape`，令 `scale[k] = τ[k] / shape`，使单位线平均到达时间等于 `τ[k]`。

### 实现要求

- 保留现有统一 `GammaUHRouter`，作为无距离数据时的简化方案；
- 新增距离感知的网格—河段坡面汇流，不在运行时追踪 flowdir；
- 网格水量必须先延迟、后汇总，避免提前求和丢失空间差异；
- 若上层水文模型已完成坡面汇流，允许直接输入 `qlat` 并跳过该过程；
- 保证单位为：距离 m、旅行时间 s、径流深率 m/s、侧向入流 m³/s。

### 验收标准

- 相同径流条件下，近河网格应早于远河网格到达河道；
- 所有网格输入水量与最终进入河道的水量守恒；
- 相同旅行时间下，新方法应退化为统一单位线汇流；
- 增加一个包含近、远两个网格的最小测试，验证到达顺序和水量守恒。
