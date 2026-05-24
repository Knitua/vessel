# STL 严格几何 QA 报告

本报告针对 `reproduction_outputs/stl_export_check/*/candidate_vasculature.stl` 做离线几何检查。

检查范围包括：水密性、边拓扑、连通分量、重复面、退化/狭长三角面、法向/体积一致性，以及 PyMeshFix 自交三角面检测。

注意：通过本 QA 仍不等于可以直接打印。制造交付前仍需先完成几何 QA 修复闭环，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。

## 总览

| 体系 | watertight | open edges | nonmanifold edges | components | duplicate face groups | self intersections | min shape quality | sliver faces <1e-4 | 结论 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `cube` | True | 0 | 0 | 1 | 0 | 0 | 0.031145164271226527 | 0 | 基础几何 QA 通过 |
| `cube_forest` | True | 0 | 0 | 1 | 0 | 0 | 0.014501929322095552 | 0 | 基础几何 QA 通过 |
| `biventricle` | True | 0 | 0 | 1 | 0 | 0 | 0.0008305848971226303 | 0 | 基础几何 QA 通过 |

## 详细结果

### `cube`

- STL: `reproduction_outputs/stl_export_check/cube/candidate_vasculature.stl`
- faces / vertices: `679368` / `339686`
- watertight: `True`; winding consistent: `True`; is_volume: `True`
- open edges: `0`; nonmanifold edges: `0`; max edge face count: `2`
- connected components: `1`; component faces: `[679368]`
- duplicate face groups: `0`
- self-intersection method: `pymeshfix.PyTMesh.select_intersecting_triangles`; intersecting triangles: `0`; elapsed: `10.421s`
- triangle area min / median / max: `1.3985278031447468e-09` / `8.76510582051406e-08` / `3.935054009688822e-07`
- edge length min / median / max: `1.378547487813685e-05` / `0.0004578999075889774` / `0.0014363178368156013`
- aspect ratio p95 / max: `1.5926747916391837` / `24.600845369597504`
- shape quality min / p01 / median: `0.031145164271226527` / `0.2623184639220905` / `0.3265029762069133`
- degenerate faces area <= 1e-16: `0`
- sliver faces quality < 1e-4 / < 1e-3: `0` / `0`

### `cube_forest`

- STL: `reproduction_outputs/stl_export_check/cube_forest/candidate_vasculature.stl`
- faces / vertices: `746196` / `373098`
- watertight: `True`; winding consistent: `True`; is_volume: `True`
- open edges: `0`; nonmanifold edges: `0`; max edge face count: `2`
- connected components: `1`; component faces: `[746196]`
- duplicate face groups: `0`
- self-intersection method: `pymeshfix.PyTMesh.select_intersecting_triangles`; intersecting triangles: `0`; elapsed: `12.83s`
- triangle area min / median / max: `1.2168954507242514e-09` / `8.224904300662881e-08` / `3.611588640582946e-07`
- edge length min / median / max: `5.4737991325932666e-05` / `0.0004426557714265975` / `0.0010609409114561963`
- aspect ratio p95 / max: `1.6086599221723419` / `6.396718648724995`
- shape quality min / p01 / median: `0.014501929322095552` / `0.2594511148209941` / `0.32664349859307207`
- degenerate faces area <= 1e-16: `0`
- sliver faces quality < 1e-4 / < 1e-3: `0` / `0`

### `biventricle`

- STL: `reproduction_outputs/stl_export_check/biventricle/candidate_vasculature.stl`
- faces / vertices: `1411274` / `705639`
- watertight: `True`; winding consistent: `True`; is_volume: `True`
- open edges: `0`; nonmanifold edges: `0`; max edge face count: `2`
- connected components: `1`; component faces: `[1411274]`
- duplicate face groups: `0`
- self-intersection method: `pymeshfix.PyTMesh.select_intersecting_triangles`; intersecting triangles: `0`; elapsed: `37.232s`
- triangle area min / median / max: `1.5282142322202078e-10` / `2.241121262407394e-07` / `9.324451314361536e-07`
- edge length min / median / max: `2.941837949574725e-06` / `0.00073334035925402` / `0.001943264937625949`
- aspect ratio p95 / max: `1.6183808726161306` / `193.01858048849365`
- shape quality min / p01 / median: `0.0008305848971226303` / `0.2592409436698991` / `0.32565573684969246`
- degenerate faces area <= 1e-16: `0`
- sliver faces quality < 1e-4 / < 1e-3: `0` / `1`
