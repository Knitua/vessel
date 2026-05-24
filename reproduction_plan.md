# svVascularize 论文软件主线复现计划

## 目标边界

本计划复现 `Rapid model-guided design of organ-scale synthetic vasculature for biomanufacturing` 在“软件与数值工作流”层面的核心链路，而不是复现湿实验、打印工艺和论文全部图表。

当前复现目标是：
- 从 domain 几何生成 vascular tree / forest。
- 导出并实际运行 0D 降阶血流模型。
- 导出 3D CFD 所需的体网格、边界面和 `fluid_simulation_*.xml`。
- 对一个小规模 tree 导出 candidate STL，并做基础几何检查。
- 审计当前 0D/3D 仿真输入输出状态，确认哪些已经求解、哪些只是下一阶段输入。

当前不作为完成标准的内容：
- 不跑论文级百万 terminal / organ-scale benchmark。
- 不跑完整 svMultiPhysics 3D CFD 求解。
- 不复现 3D 打印、OCT、PIV、LIVE/DEAD、7 天灌流实验。
- 不声明 STL 已经是可直接打印的制造交付件。

## 当前环境

当前复现使用 `vessel` 环境：
- Python: `3.10.20`
- PyVista: `0.48.4`
- MMG: `svv/utils/remeshing/Linux/x86_64/{mmg2d_O3,mmg3d_O3,mmgs_O3}`
- 0D solver: `svv/utils/solvers/0D/Linux/x86_64/svzerodsolver`
- svMultiPhysics: 当前 `PATH` 中未找到，因此 3D CFD 尚未实际启动。
- headless 设置：`SVV_TELEMETRY_DISABLED=1`、`PYVISTA_OFF_SCREEN=true`、`QT_QPA_PLATFORM=offscreen`

复现入口脚本：
- [`run_reproduction_plan.py`](/root/autodl-tmp/svVascularize/run_reproduction_plan.py)

输出目录：
- [`reproduction_outputs/`](/root/autodl-tmp/svVascularize/reproduction_outputs)
- 汇总文件：[`reproduction_outputs/summary.json`](/root/autodl-tmp/svVascularize/reproduction_outputs/summary.json)
- 中文输出说明：[`reproduction_outputs_README.md`](/root/autodl-tmp/svVascularize/reproduction_outputs_README.md)

## 已完成的四个目标体系

| 体系 | 目标 | 当前状态 |
|---|---|---|
| `cube` | 小立方体 domain 中生成单树，跑 0D，导出 3D 输入 | 已完成 |
| `cube_3d` | 更小的 cube 单树，验证 0D + 3D 输入导出 | 已完成 |
| `cube_forest` | 两棵树组成 connected forest，跑 0D，导出 3D 输入 | 已完成 |
| `biventricle` | 复杂 STL domain 中生成 tree，跑 0D，导出 3D 输入 | 已完成 |

规模和输出摘要：

| 体系 | 几何/网络规模 | 0D 状态 | 3D 输入导出状态 | 3D 体网格 cells | 3D surface cells |
|---|---:|---|---|---:|---:|
| `cube` | 单树，7 segments，4 terminals | 成功 | 成功 | 2,805,147 | 679,368 |
| `cube_3d` | 单树，3 segments，2 terminals | 成功 | 成功 | 1,758,519 | 411,664 |
| `cube_forest` | 1 network，2 棵树，每棵 3 segments，1 个 connection | 成功 | 成功 | 3,016,130 | 746,196 |
| `biventricle` | 复杂 STL domain，单树，3 segments，2 terminals | 成功 | 成功 | 5,755,838 | 1,411,274 |

## 当前执行链路

`run_reproduction_plan.py` 按以下顺序执行：

1. 建立 cube domain：`Domain(pv.Cube().triangulate()) -> create() -> solve() -> build(resolution=15)`。
2. 在 cube domain 中生成 `cube` tree：`set_domain() -> set_root() -> n_add(3)`。
3. 导出 `cube` centerlines、0D 输入、0D 结果、3D mesh/XML。
4. 从 `cube` tree 直接导出 candidate STL，并做水密性、连通性、最小半径/管径检查。
5. 在同一 cube domain 中生成 `cube_3d` 小树：`n_add(1)`。
6. 导出 `cube_3d` centerlines、0D 输入、0D 结果、3D mesh/XML。
7. 在 cube domain 中生成 `cube_forest`：`Forest(..., n_trees_per_network=[2]) -> set_roots() -> add(1) -> connect()`。
8. 导出并运行 `cube_forest` 的 0D，导出 connected forest 的 3D mesh/XML。
9. 读取 `biventricle_tissue.stl`，建立非凸 domain：`create() -> solve() -> build(resolution=12)`。
10. 在 biventricle domain 中生成小树：`n_add(1)`。
11. 导出 `biventricle` centerlines、0D 输入、0D 结果、3D mesh/XML。

## 输出文件含义

| 文件类型 | 含义 | 当前用途 |
|---|---|---|
| `solver_0d.in` | 0D solver 输入 | 用于 `svZeroDSolver` |
| `output.csv` | 0D solver 输出结果 | 已经是真正的 0D 计算结果 |
| `geom.csv` | 0D/几何辅助数据 | 辅助检查和后处理 |
| `centerlines.vtp` | 中心线几何 | 可视化、后处理、几何核查 |
| `fluid_msh_*.vtu` | 3D volume mesh | 下一阶段 3D CFD 输入 |
| `mesh-surfaces/*.vtp` | wall/cap/lumen 边界面 | 下一阶段 3D CFD 边界条件输入 |
| `fluid_simulation_*.xml` | svMultiPhysics solver 配置 | 下一阶段 3D CFD 输入 |
| `candidate_vasculature.stl` | 候选制造几何 | 几何质量检查，不等于可直接打印 |

## 0D、3D 和 STL 的关系

`0D` 是降阶网络模型。它把血管段近似为 lumped-parameter 元件，用来快速计算压力和流量。当前 `output.csv` 是已经跑出来的 0D 数值结果。

`3D` 当前指 3D CFD 输入导出，包括 volume mesh、boundary surfaces 和 XML。它还不是 3D CFD 结果；真正的 3D 结果需要下一阶段用 svMultiPhysics 求解。

`STL` 是制造/几何方向的候选表面文件。当前 STL 覆盖 `cube`、`cube_forest` 和 `biventricle`：`cube`/`biventricle` 从 tree 的 `centerline + radius` 直接导出 watertight surface，`cube_forest` 从 connected forest solid surface 导出。它们不是从 3D CFD 结果转出来的，也没有包含仿真验证。

面向制造交付的推荐顺序是：先做严格几何 QA，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。当前 STL 不能直接作为 print-ready model。

## 验收结果

当前已通过的验收：
- `cube` 同时有 0D solver 输出和 3D mesh/XML。
- `cube_3d` 同时有 0D solver 输出和 3D mesh/XML。
- `cube_forest` 完成 connected forest、0D solver 输出和 connected 3D mesh/XML。
- `biventricle` 完成复杂 STL domain 中的小树生成、0D solver 输出和 3D mesh/XML。
- `stl_export_check` 为 `cube`、`cube_forest`、`biventricle` 生成 candidate STL，且基础检查显示 watertight、1 个 connected component、open edges 为 0。
- `stl_qa` 已完成三个 STL 的严格几何 QA：0 open edges、0 nonmanifold edges、0 duplicate face groups、0 self-intersecting triangles；`biventricle` 仍有 1 个 `shape_quality < 1e-3` 的狭长三角面，制造前应局部复核。
- `simulation_audit` 已完成仿真侧审计：0D 四个体系均有真实 `output.csv`；3D 四个体系均有 mesh/XML 且所有 face 已有边界条件，但没有 svMultiPhysics 可执行程序和 3D `result*` 输出。
- 相关测试子集通过：`27 passed`。

测试命令：

```bash
SVV_TELEMETRY_DISABLED=1 PYVISTA_OFF_SCREEN=true QT_QPA_PLATFORM=offscreen \
/root/miniconda3/envs/vessel/bin/python -m pytest -q \
  test/test_forest_assignment.py \
  test/test_dmn_io.py \
  test/test_voronoi_sampling.py \
  test/test_solver_0d_selector.py \
  test/test_zerod_export_solver_opt_in.py \
  test/test_spline_export_files.py \
  test/test_remesh_surface.py
```

复现脚本运行命令：

```bash
SVV_TELEMETRY_DISABLED=1 PYVISTA_OFF_SCREEN=true QT_QPA_PLATFORM=offscreen \
/root/miniconda3/envs/vessel/bin/python run_reproduction_plan.py
```

复用已有重型 3D 输出时：

```bash
SVV_REPRO_REUSE_OUTPUTS=1 SVV_TELEMETRY_DISABLED=1 PYVISTA_OFF_SCREEN=true QT_QPA_PLATFORM=offscreen \
/root/miniconda3/envs/vessel/bin/python run_reproduction_plan.py
```

## 和论文算法的一致性说明

当前复现使用的是仓库自身的公开接口和算法流程，包括 `Domain`、`Tree`、`Forest`、`TreeConnection`、0D export、3D export、MMG remeshing、TetGen meshing 等。

但当前不是严格复现论文的全规模结果：
- 规模显著缩小，只做小树和小 forest。
- 没有复现百万 terminal、批量 200+ geometry、论文计时 benchmark。
- 没有跑完整 3D CFD，也没有做 0D-3D coupling。
- 没有做打印和湿实验验证。

为了让当前仓库在本机环境跑通，我做了若干兼容性修复和小规模鲁棒性兜底。具体见：
- [`reproduction_fix_report.md`](/root/autodl-tmp/svVascularize/reproduction_fix_report.md)

## 下一阶段建议

若目标是仿真：
1. 从 `reproduction_outputs/simulation_audit/simulation_audit.md` 复核当前 0D/3D 审计结果。
2. 安装并验证 svMultiPhysics。
3. 检查 `fluid_simulation_*.xml` 中的材料参数、时间步、边界条件和输出设置。
4. 先从 `cube_3d` 跑短时间 3D CFD smoke solve。
5. 再跑 `cube`、`cube_forest`、`biventricle`。
6. 对比 0D 和 3D 的入口/出口流量、压力范围和守恒误差。

若目标是制造/打印：
1. 从 `reproduction_outputs/stl_qa/strict_geometry_qa.md` 和 `reproduction_outputs/stl_qa/strict_geometry_qa.json` 复核当前几何 QA 结果。
2. 对 `biventricle` 的低质量狭长三角面做局部网格/曲率复核。
3. 做 0D/3D 仿真验证，确认灌流、压降、流量分布和局部流场风险。
4. 根据具体打印工艺加入最小管径、支撑、清洗、壁厚、出口连通性、材料收缩和后处理约束。
5. 不建议直接把当前 STL 当作最终打印文件。
