# svVascularize 复现修复说明

本文档记录本轮为了跑通 `reproduction_plan.md` 做过的主要修复，以及这些修复和论文算法复现之间的关系。

## 结论先行

当前结果不是“严格复现论文全部规模和全部实验”，而是“使用仓库原始算法接口，跑通小规模软件主线”。

已经跑通：
- `cube`: 0D 求解 + 3D mesh/XML 导出。
- `cube_3d`: 0D 求解 + 3D mesh/XML 导出。
- `cube_forest`: connected forest + 0D 求解 + 3D mesh/XML 导出。
- `biventricle`: 复杂 STL domain 中生成 tree + 0D 求解 + 3D mesh/XML 导出。
- `stl_export_check`: 从小 `cube` tree 直接导出 candidate STL，并做基础几何检查。
- `stl_qa`: 对 `cube`、`cube_forest`、`biventricle` 三个 candidate STL 完成严格几何 QA，报告见 `reproduction_outputs/stl_qa/strict_geometry_qa.md`。
- `simulation_audit`: 对四个体系完成 0D/3D 仿真侧审计，报告见 `reproduction_outputs/simulation_audit/simulation_audit.md`。

没有完成：
- 没有跑论文级百万 terminal / organ-scale benchmark。
- 没有跑 svMultiPhysics 完整 3D CFD 求解。
- 当前机器 `PATH` 中未找到 svMultiPhysics，可用状态尚未验证。
- 没有复现 3D 打印、OCT、PIV、LIVE/DEAD、灌流实验。
- 当前 STL 不是可直接打印的最终制造文件；面向制造交付时，应先做严格几何 QA，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。

## 问题类型概览

| 类型 | 说明 | 对算法复现的影响 |
|---|---|---|
| 环境/依赖问题 | 当前机器缺少或不兼容部分依赖、二进制、headless GUI 设置 | 需要补齐环境后才能运行 |
| PyVista API 兼容问题 | 新版 PyVista 的 cell quality API 和旧代码预期不一致 | 不改变算法，只是兼容新版库 |
| 输出路径问题 | `Tree.export_solid()` 没有保证输出目录存在 | 不改变算法，只是让导出稳定 |
| mesh 属性保存问题 | 旧代码依赖 `mesh.hsize` 动态属性，但 PyVista 对象序列化/过滤后属性不可靠 | 不改变算法，改为用 `cell_data["hsize"]` 传递 |
| forest sparse assignment 问题 | 小规模 forest terminal 数很少时 sparse matching 可能没有 full matching | 增加 dense assignment 兜底，属于鲁棒性修复 |
| connected forest 3D solid 导出问题 | `TreeConnection.export_solid()` 某些分支返回 `None` 或写错内部列表 | 属于真实代码缺陷修复 |
| 3D face 映射问题 | connected forest 体网格后，用原 surface 映射 volume face 会出现 ID 不一致 | 改用 volume mesh 提取 surface，属于导出链路修复 |
| 3D wall/lumen 边界条件缺口 | `lumen_*` 血管壁面已导出为 face，但旧 XML 没有为其添加 no-slip BC | 属于 3D 输入完整性修复 |

## 主要修复清单

### 1. PyVista cell quality 兼容层

新增文件：
- [`svv/utils/pyvista_compat.py`](/root/autodl-tmp/svVascularize/svv/utils/pyvista_compat.py)

问题：
- 原代码多处调用 `mesh.compute_cell_quality(...).cell_data["CellQuality"]`。
- 当前 PyVista 版本中相关 API/输出字段存在兼容差异，可能导致质量检查取不到预期字段。

修复：
- 增加 `compute_cell_quality()` 兼容函数。
- 如果新版 API 输出字段名不是 `CellQuality`，则补齐 legacy-compatible `CellQuality` cell array。

影响：
- 不改变 remeshing、surface extraction 或 mesh quality 算法。
- 只是让原有质量检查逻辑能在当前 PyVista 版本下运行。

涉及文件：
- [`svv/tree/export/export_solid.py`](/root/autodl-tmp/svVascularize/svv/tree/export/export_solid.py)
- [`svv/simulation/utils/extract_faces.py`](/root/autodl-tmp/svVascularize/svv/simulation/utils/extract_faces.py)
- [`svv/simulation/utils/boundary_layer.py`](/root/autodl-tmp/svVascularize/svv/simulation/utils/boundary_layer.py)
- [`svv/utils/remeshing/remesh.py`](/root/autodl-tmp/svVascularize/svv/utils/remeshing/remesh.py)

### 2. `Domain.get_interior_points(method="voronoi")` 补齐

涉及文件：
- [`svv/domain/domain.py`](/root/autodl-tmp/svVascularize/svv/domain/domain.py)

问题：
- 代码路径中存在 `method="voronoi"` 的使用预期，但原函数没有明确处理这个分支。
- 在某些采样/测试路径下会退到不合适的默认逻辑。

修复：
- 增加 `voronoi` 分支，使用 domain mesh cell centers 作为候选点。
- 支持 `implicit_range`、tree 距离 threshold、`volume_threshold` 过滤。

影响：
- 这是补齐已有接口预期，不是新增论文外算法目标。
- 当前复现规模小，主要用于让采样逻辑和相关测试稳定。

### 3. `Tree.export_solid()` 输出目录修复

涉及文件：
- [`svv/tree/tree.py`](/root/autodl-tmp/svVascularize/svv/tree/tree.py)

问题：
- `export_solid()` 写文件时没有保证 `outdir` 存在，某些导出路径会失败。

修复：
- 在导出前执行 `os.makedirs(outdir, exist_ok=True)`。

影响：
- 不改变几何生成算法，只修复文件系统导出行为。

### 4. watertight solid 的 `hsize` 保存位置修复

涉及文件：
- [`svv/tree/export/export_solid.py`](/root/autodl-tmp/svVascularize/svv/tree/export/export_solid.py)

问题：
- 原代码在 poor-quality element 修复前后对 `hsize` 的保存时机不稳。
- mesh repair 后动态属性可能丢失。

修复：
- 在 repair 前先从 `model.cell_data["hsize"][0]` 保存 `hsize`。
- repair 后重新写回 `cell_data["hsize"]`。

影响：
- 不改变 vessel solid 生成方法，只保证后续 remeshing/3D export 能拿到目标网格尺度。

### 5. forest terminal assignment 兜底

涉及文件：
- [`svv/forest/connect/assign.py`](/root/autodl-tmp/svVascularize/svv/forest/connect/assign.py)
- [`test/test_forest_assignment.py`](/root/autodl-tmp/svVascularize/test/test_forest_assignment.py)

问题：
- `cube_forest` 小规模情况下，每棵树 terminal 数很少。
- 默认 sparse candidate graph 可能不足以形成 full bipartite matching。
- 原逻辑遇到 sparse matching 失败后直接返回 `None, None`，导致 `Forest.connect()` 后网络拓扑不完整。
- 结果是 0D solver 会遇到 singular system。

修复：
- sparse matching 失败时，退回 dense terminal cost matrix。
- 对 convex domain 使用 Euclidean distance。
- 对非凸 domain 先用已有 geodesic cost；缺失项再用 Euclidean distance 兜底。
- 如果 dense assignment 选中的 pair 不在原 sparse path map 里，就补一个线性 spline connection function。

影响：
- 这是小规模鲁棒性修复，不是论文中大规模 optimal assignment 的性能复现。
- 对当前 `cube_forest`，它保证 connected forest 拓扑完整，从而 0D 和 3D 导出可以继续。
- 严格大规模复现时，应调大 neighbor 搜索或使用论文级配置，而不是依赖 dense fallback 作为性能路径。

回归测试：
- 新增测试强制 sparse matching 抛错，并验证 fallback 会补齐 connection functions。

### 6. connected forest solid export 修复

涉及文件：
- [`svv/forest/connect/tree_connection.py`](/root/autodl-tmp/svVascularize/svv/forest/connect/tree_connection.py)

问题 1：
- `p.shape[1] == 3` 分支里把 spline 追加到了错误的列表，导致后续半径/路径插值状态不一致。

修复 1：
- 从 `interp_xyz.append(...)` 改为 `network_xyz.append(...)`。

问题 2：
- 当 `junction_smoothing=False` 时，原代码直接 `network_solids.append(None)`。
- 这导致 connected forest 3D export 拿不到 surface mesh，后续保存或 tetrahedralize 会失败。

修复 2：
- 对生成的 tubes 做 `union_tubes(...)`，失败时用 `pv.merge(tubes)` 兜底。
- 单 tube 情况直接使用该 tube surface。
- 写入 `cell_data["hsize"]`。
- 最终 append 实际 `model`，而不是 `None`。

影响：
- 这是 connected forest 3D 导出链路的真实缺陷修复。
- `pv.merge(tubes)` fallback 是几何导出的兜底，适用于当前小规模 connected forest；严格高质量复杂 junction 仍应优先使用 union/smoothing 路线并做几何 QA。

### 7. connected forest 3D simulation 分支修复

涉及文件：
- [`svv/simulation/simulation.py`](/root/autodl-tmp/svVascularize/svv/simulation/simulation.py)

问题 1：
- 原代码直接读 `fluid_surface.hsize` 动态属性。
- PyVista filter/tetgen/repair 后该属性可能不存在。

修复 1：
- 优先读 `getattr(fluid_surface, "hsize", None)`。
- 如果没有，则读 `fluid_surface.cell_data["hsize"][0]`。
- 仍没有时，用 `fluid_surface.length / 100.0` 作为保底尺度。

问题 2：
- 原代码在 connected forest 体网格生成后，仍可能用原始 surface 做 face extraction。
- TetGen 后 volume mesh 的点/面 ID 和原始 surface 不完全一致，导致 `extract_faces()` 映射 `GlobalElementID` 时出现 KeyError。

修复 2：
- 体网格生成成功后，统一使用 `fluid_volume.extract_surface()` 作为后续 face extraction surface。
- 这样边界面与 volume mesh cell id 对应关系一致。

问题 3：
- boundary layer / wall layer 逻辑在没有可用 wall/lumen 时可能继续访问空列表。

修复 3：
- 增加空 wall/lumen 检查；没有可用 wall 时跳过相应 layer 并记录 `None`。

影响：
- 不改变 3D CFD 求解算法，因为当前阶段只导出 mesh/XML，没有运行 svMultiPhysics。
- 这是让 connected forest 的 3D 输入导出在当前 PyVista/TetGen 组合下可完成。

## 复现脚本说明

新增脚本：
- [`run_reproduction_plan.py`](/root/autodl-tmp/svVascularize/run_reproduction_plan.py)

脚本做了什么：
- 固定随机种子。
- 构建 `cube`、`cube_3d`、`cube_forest`、`biventricle` 四个小规模体系。
- 对每个主体系导出并运行 0D。
- 对每个主体系导出 3D mesh/XML。
- 对 `cube`、`cube_forest`、`biventricle` 额外导出 candidate STL。
- 把文件路径、网格规模、状态写入 `reproduction_outputs/summary.json`。

脚本的简化：
- 使用非常小的 `n_add`，不是论文级规模。
- 3D 只导出输入，不跑 svMultiPhysics。
- 0D 只跑短周期、少时间点，用于验证链路。
- STL 只作为 candidate manufacturing geometry 导出，不代表已经完成 organ-scale 打印模型或制造工艺验证。后续应先做严格几何 QA，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。

## 是否严格复现文章算法

答案需要分层看。

严格使用的部分：
- 使用仓库原生 `Domain`、`Tree`、`Forest`、`Simulation` 接口。
- 使用仓库已有 tree growth、forest connection、solid export、0D export、3D export 流程。
- 使用仓库的 MMG/TetGen 网格链路。
- 使用仓库提供/配置的 `svZeroDSolver` 实际求解 0D。

不严格的部分：
- 没有使用论文中的大规模参数。
- 没有复现论文 reported runtime benchmark。
- 没有跑完整 3D CFD 或 0D-3D coupling。
- `cube_forest` 的 sparse assignment 失败时加入 dense fallback，这是为了小规模鲁棒性，不是论文大规模性能路径。
- connected forest 3D export 中 `pv.merge(tubes)` 是 union 失败时的兜底，属于工程化补救，需要后续更严格几何 QA。

因此当前应表述为：

> 使用论文代码仓的核心接口与算法流程，完成了小规模软件主线复现和输入/输出验证；为适配当前环境与修复仓库运行缺陷，加入了若干兼容性修复和小规模兜底。它不是论文全规模、全物理、全实验结果的严格复现。

## 当前可引用的结果文件

输出说明：
- [`reproduction_outputs_README.md`](/root/autodl-tmp/svVascularize/reproduction_outputs_README.md)

原始机器可读摘要：
- [`reproduction_outputs/summary.json`](/root/autodl-tmp/svVascularize/reproduction_outputs/summary.json)

关键输出：
- `reproduction_outputs/cube/tree_0d/output.csv`
- `reproduction_outputs/cube/tree_3d/fluid_simulation_0-0.xml`
- `reproduction_outputs/cube_3d/tree_0d/output.csv`
- `reproduction_outputs/cube_3d/tree_3d/fluid_simulation_0-0.xml`
- `reproduction_outputs/cube_forest/forest_0d/output.csv`
- `reproduction_outputs/cube_forest/forest_3d/fluid_simulation_0-0.xml`
- `reproduction_outputs/biventricle/tree_0d/output.csv`
- `reproduction_outputs/biventricle/tree_3d/fluid_simulation_0-0.xml`
- `reproduction_outputs/stl_export_check/cube/candidate_vasculature.stl`
- `reproduction_outputs/stl_export_check/cube_forest/candidate_vasculature.stl`
- `reproduction_outputs/stl_export_check/biventricle/candidate_vasculature.stl`
- `reproduction_outputs/stl_qa/strict_geometry_qa.md`
- `reproduction_outputs/stl_qa/strict_geometry_qa.json`
- `reproduction_outputs/simulation_audit/simulation_audit.md`
- `reproduction_outputs/simulation_audit/simulation_audit.json`
