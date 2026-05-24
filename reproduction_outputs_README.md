# 复现输出说明

`run_reproduction_plan.py` 会把小规模软件复现实验的结果写入 `reproduction_outputs/`。

当前结果用于验证 `svVascularize` 的核心软件流程：
- `cube`: `Domain -> Tree -> 0D solver -> 3D mesh/XML -> centerlines`
- `cube_3d`: `Tree -> 0D solver -> 3D mesh/XML`
- `cube_forest`: `Forest -> connect -> 0D solver -> 3D mesh/XML`
- `biventricle`: `STL Domain -> Tree -> 0D solver -> 3D mesh/XML -> centerlines`
- `stl_export_check`: `cube`、`cube_forest`、`biventricle` 的候选 STL 导出和几何检查

## 当前阶段结论

现在已经完成的是小规模 workflow validation，不是论文中的百万 outlet / organ-scale benchmark。

- `0D`: 已真正调用 `svZeroDSolver` 跑完，并生成 `output.csv`。这是降阶血流网络模型的数值计算结果。
- `3D`: 已生成体网格、边界面和 `fluid_simulation_*.xml`。这些是 svMultiPhysics 的 3D CFD 输入文件；已补齐 `lumen_*` 血管壁面的 no-slip 边界条件，但本机尚未发现 svMultiPhysics 可执行程序，也没有 3D `result*` 输出文件，所以还没有完成 3D CFD 求解。
- `STL`: 已生成 `cube`、`cube_forest`、`biventricle` 三个候选制造几何并检查水密性；它们还不是经过制造约束验证的可直接打印交付件。

## 四个体系规模和状态

| 体系 | 几何/网络规模 | 0D 状态 | 3D 输入导出状态 | 3D 体网格 cells | 3D surface cells |
|---|---:|---|---|---:|---:|
| `cube` | 单树，7 segments，4 terminals | 成功 | 成功 | 2,805,147 | 679,368 |
| `cube_3d` | 单树，3 segments，2 terminals | 成功 | 成功 | 1,758,519 | 411,664 |
| `cube_forest` | 1 network，2 棵树，每棵 3 segments，1 个 connection | 成功 | 成功 | 3,016,130 | 746,196 |
| `biventricle` | 复杂 STL domain，单树，3 segments，2 terminals | 成功 | 成功 | 5,755,838 | 1,411,274 |

## Domain 和中心线规模

| 体系 | domain points | patches | boundary cells | domain mesh cells | convexity | centerline points | centerline cells |
|---|---:|---:|---:|---:|---:|---:|---:|
| `cube` | 24 | 24 | 12 | 6 | 1.0 | 400 | 399 |
| `cube_3d` | 复用 cube domain | 复用 cube domain | 复用 cube domain | 复用 cube domain | 复用 cube domain | 200 | 199 |
| `cube_forest` | 复用 cube domain | 复用 cube domain | 复用 cube domain | 复用 cube domain | 复用 cube domain | 未单独导出 | 未单独导出 |
| `biventricle` | 19,368 | 19,364 | 37,014 | 79,787 | 0.6764824779 | 200 | 199 |

## 输出文件用途

| 文件类型 | 示例 | 用途 | 是否是仿真输入 | 是否是仿真结果 |
|---|---|---|---|---|
| `solver_0d.in` | `tree_0d/solver_0d.in` | 0D solver 输入 | 是，0D | 否 |
| `output.csv` | `tree_0d/output.csv` | 0D solver 输出结果 | 否 | 是，0D |
| `geom.csv` | `tree_0d/geom.csv` | 0D/几何辅助信息 | 部分是 | 否 |
| `centerlines.vtp` | `centerlines/centerlines.vtp` | 中心线几何、可视化、后处理 | 可作为辅助输入 | 否 |
| `fluid_msh_*.vtu` | `tree_3d/mesh/.../*.vtu` | 3D volume mesh | 是，3D | 否 |
| `mesh-surfaces/*.vtp` | `mesh-surfaces/wall_0.vtp` | 3D 边界面，用于边界条件 | 是，3D | 否 |
| `fluid_simulation_*.xml` | `tree_3d/fluid_simulation_0-0.xml` | svMultiPhysics solver 配置 | 是，3D | 否 |
| `candidate_vasculature.stl` | `stl_export_check/*/candidate_vasculature.stl` | 候选制造/几何交付 | 否 | 否 |

## 关键输出路径

| 体系 | 0D 输入 | 0D 结果 | 3D XML | 3D mesh 目录 |
|---|---|---|---|---|
| `cube` | `reproduction_outputs/cube/tree_0d/solver_0d.in` | `reproduction_outputs/cube/tree_0d/output.csv` | `reproduction_outputs/cube/tree_3d/fluid_simulation_0-0.xml` | `reproduction_outputs/cube/tree_3d/mesh` |
| `cube_3d` | `reproduction_outputs/cube_3d/tree_0d/solver_0d.in` | `reproduction_outputs/cube_3d/tree_0d/output.csv` | `reproduction_outputs/cube_3d/tree_3d/fluid_simulation_0-0.xml` | `reproduction_outputs/cube_3d/tree_3d/mesh` |
| `cube_forest` | `reproduction_outputs/cube_forest/forest_0d/solver_0d.in` | `reproduction_outputs/cube_forest/forest_0d/output.csv` | `reproduction_outputs/cube_forest/forest_3d/fluid_simulation_0-0.xml` | `reproduction_outputs/cube_forest/forest_3d/mesh` |
| `biventricle` | `reproduction_outputs/biventricle/tree_0d/solver_0d.in` | `reproduction_outputs/biventricle/tree_0d/output.csv` | `reproduction_outputs/biventricle/tree_3d/fluid_simulation_0-0.xml` | `reproduction_outputs/biventricle/tree_3d/mesh` |

## STL 候选几何检查

| 体系 | STL 路径 | surface cells | surface points | open edges | PyVista 水密性 | Trimesh 水密性 | connected components | min radius | min diameter |
|---|---|---:|---:|---:|---|---|---:|---:|---:|
| `cube` | `reproduction_outputs/stl_export_check/cube/candidate_vasculature.stl` | 679,368 | 339,686 | 0 | true | true | 1 | 0.0027793745353650564 | 0.005558749070730113 |
| `cube_forest` | `reproduction_outputs/stl_export_check/cube_forest/candidate_vasculature.stl` | 746,196 | 373,098 | 0 | true | true | 1 | 0.0033527637561611384 | 0.006705527512322277 |
| `biventricle` | `reproduction_outputs/stl_export_check/biventricle/candidate_vasculature.stl` | 1,411,274 | 705,639 | 0 | true | true | 1 | 0.0044679034197906975 | 0.008935806839581395 |

以上 STL 均为 candidate manufacturing geometry。旧的 STL 导出摘要中 `self_intersection_check` 字段为 `not_available`，不能只依据该摘要声明已经完成严格自交检测，也不能说明可以直接开始打印。

已补充严格几何 QA：
- Markdown 报告：`reproduction_outputs/stl_qa/strict_geometry_qa.md`
- JSON 原始结果：`reproduction_outputs/stl_qa/strict_geometry_qa.json`

QA 结果显示三个 STL 均为 watertight、0 open edges、0 nonmanifold edges、1 个 connected component、0 duplicate face groups、PyMeshFix 自交三角面检测为 0。剩余注意点：`biventricle` 有 1 个 `shape_quality < 1e-3` 的狭长三角面，制造前仍建议按具体打印工艺做局部网格/曲率复核。

面向制造交付的推荐顺序是：先做严格几何 QA，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。

## 和下一阶段的关系

当前 3D 输出已经是 svMultiPhysics 可读取方向的输入产物，但还没有完成 3D CFD 求解。下一阶段如果要进入真正仿真，需要：

1. 检查 `fluid_simulation_*.xml` 中的材料参数、时间步、边界条件和输出设置是否符合目标实验。
2. 用 svMultiPhysics 读取 XML 和 mesh 运行 3D CFD。
3. 检查 svMultiPhysics 的收敛、流量守恒、压力范围和壁面边界标记。
4. 若目标是制造交付，则先做严格几何 QA，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。

完整原始摘要见 `reproduction_outputs/summary.json`，里面包含精确文件路径、mesh cell counts、运行时间和状态字段。

## 仿真审计

已补充仿真侧审计报告：
- Markdown 报告：`reproduction_outputs/simulation_audit/simulation_audit.md`
- JSON 原始结果：`reproduction_outputs/simulation_audit/simulation_audit.json`

审计结果显示：
- 0D 四个体系均已有真实 `output.csv`，最大相对流量守恒误差约 `4.557e-08`。
- 3D 四个体系的 mesh/XML 文件存在，所有导出的 face 都已有边界条件。
- 当前机器 `PATH` 中未找到 `svMultiPhysics`，也未发现 3D `result*` 输出文件；因此 3D 仍处于输入准备阶段。
