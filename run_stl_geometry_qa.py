"""Run geometry QA checks for reproduction STL candidates."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh
from pymeshfix import _meshfix


REPO_ROOT = Path(__file__).resolve().parent
OUT_DIR = REPO_ROOT / "reproduction_outputs" / "stl_qa"
STL_CASES = {
    "cube": REPO_ROOT / "reproduction_outputs" / "stl_export_check" / "cube" / "candidate_vasculature.stl",
    "cube_forest": REPO_ROOT / "reproduction_outputs" / "stl_export_check" / "cube_forest" / "candidate_vasculature.stl",
    "biventricle": REPO_ROOT / "reproduction_outputs" / "stl_export_check" / "biventricle" / "candidate_vasculature.stl",
}


def _safe_float(value):
    value = float(value)
    if math.isfinite(value):
        return value
    return None


def _summary_stats(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"min": None, "p01": None, "p05": None, "median": None, "p95": None, "max": None}
    return {
        "min": _safe_float(np.min(values)),
        "p01": _safe_float(np.percentile(values, 1)),
        "p05": _safe_float(np.percentile(values, 5)),
        "median": _safe_float(np.percentile(values, 50)),
        "p95": _safe_float(np.percentile(values, 95)),
        "max": _safe_float(np.max(values)),
    }


def _edge_counts(faces: np.ndarray) -> dict:
    edges = np.vstack(
        [
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ]
    )
    edges.sort(axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return {
        "unique_edges": int(unique.shape[0]),
        "boundary_edges": int(np.sum(counts == 1)),
        "manifold_edges": int(np.sum(counts == 2)),
        "nonmanifold_edges": int(np.sum(counts > 2)),
        "max_edge_face_count": int(np.max(counts)) if counts.size else 0,
    }


def _triangle_quality(points: np.ndarray, faces: np.ndarray) -> dict:
    tri = points[faces]
    edge_lengths = np.stack(
        [
            np.linalg.norm(tri[:, 1] - tri[:, 0], axis=1),
            np.linalg.norm(tri[:, 2] - tri[:, 1], axis=1),
            np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1),
        ],
        axis=1,
    )
    area = trimesh.triangles.area(tri)
    perimeter = np.sum(edge_lengths, axis=1)
    longest = np.max(edge_lengths, axis=1)
    shortest = np.min(edge_lengths, axis=1)
    aspect_ratio = np.divide(longest, shortest, out=np.full_like(longest, np.inf), where=shortest > 0)
    # Normalized shape quality: 1 for equilateral, approaches 0 for degenerate.
    shape_quality = np.divide(
        4.0 * np.sqrt(3.0) * area,
        perimeter * perimeter,
        out=np.zeros_like(area),
        where=perimeter > 0,
    )
    return {
        "area": _summary_stats(area),
        "edge_length": _summary_stats(edge_lengths.reshape(-1)),
        "aspect_ratio": _summary_stats(aspect_ratio),
        "shape_quality": _summary_stats(shape_quality),
        "degenerate_faces_area_le_1e-16": int(np.sum(area <= 1e-16)),
        "sliver_faces_quality_lt_1e-4": int(np.sum(shape_quality < 1e-4)),
        "sliver_faces_quality_lt_1e-3": int(np.sum(shape_quality < 1e-3)),
    }


def _duplicate_faces(faces: np.ndarray) -> int:
    canonical = np.sort(faces, axis=1)
    _, counts = np.unique(canonical, axis=0, return_counts=True)
    return int(np.sum(counts > 1))


def _self_intersections(points: np.ndarray, faces: np.ndarray) -> dict:
    start = time.perf_counter()
    tmesh = _meshfix.PyTMesh()
    tmesh.load_array(points.astype(np.float64, copy=False), faces.astype(np.int32, copy=False))
    intersecting = tmesh.select_intersecting_triangles()
    elapsed = time.perf_counter() - start
    return {
        "status": "completed",
        "method": "pymeshfix.PyTMesh.select_intersecting_triangles",
        "intersecting_triangle_count": int(len(intersecting)),
        "has_self_intersections": bool(len(intersecting) > 0),
        "elapsed_sec": round(elapsed, 3),
    }


def qa_one(name: str, path: Path) -> dict:
    start = time.perf_counter()
    mesh = trimesh.load_mesh(str(path), force="mesh", process=True)
    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError(f"{path} did not load as a Trimesh")
    mesh.merge_vertices()
    mesh.update_faces(mesh.unique_faces())
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()

    points = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    components = mesh.split(only_watertight=False)
    pv_mesh = pv.read(path)

    result = {
        "name": name,
        "path": str(path.relative_to(REPO_ROOT)),
        "file_size_bytes": int(path.stat().st_size),
        "vertices": int(points.shape[0]),
        "faces": int(faces.shape[0]),
        "trimesh_is_watertight": bool(mesh.is_watertight),
        "trimesh_is_winding_consistent": bool(mesh.is_winding_consistent),
        "trimesh_euler_number": int(mesh.euler_number),
        "pyvista_open_edges": int(pv_mesh.n_open_edges),
        "connected_components": int(len(components)),
        "component_faces": [int(component.faces.shape[0]) for component in components],
        "bounds": [[_safe_float(v) for v in row] for row in np.asarray(mesh.bounds)],
        "extents": [_safe_float(v) for v in np.asarray(mesh.extents)],
        "surface_area": _safe_float(mesh.area),
        "volume": _safe_float(mesh.volume),
        "is_volume": bool(mesh.is_volume),
        "edge_counts": _edge_counts(faces),
        "duplicate_face_groups": _duplicate_faces(faces),
        "triangle_quality": _triangle_quality(points, faces),
    }
    result["self_intersections"] = _self_intersections(points, faces)
    result["elapsed_sec"] = round(time.perf_counter() - start, 3)
    return result


def _format_bool(value: bool) -> str:
    return "通过" if value else "未通过"


def write_markdown(results: dict, path: Path) -> None:
    lines = [
        "# STL 严格几何 QA 报告",
        "",
        "本报告针对 `reproduction_outputs/stl_export_check/*/candidate_vasculature.stl` 做离线几何检查。",
        "",
        "检查范围包括：水密性、边拓扑、连通分量、重复面、退化/狭长三角面、法向/体积一致性，以及 PyMeshFix 自交三角面检测。",
        "",
        "注意：通过本 QA 仍不等于可以直接打印。制造交付前仍需先完成几何 QA 修复闭环，再做 0D/3D 仿真验证，最后结合具体打印工艺做制造约束检查。",
        "",
        "## 总览",
        "",
        "| 体系 | watertight | open edges | nonmanifold edges | components | duplicate face groups | self intersections | min shape quality | sliver faces <1e-4 | 结论 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]

    for name, result in results["cases"].items():
        edge_counts = result["edge_counts"]
        quality = result["triangle_quality"]
        intersections = result["self_intersections"]
        passed = (
            result["trimesh_is_watertight"]
            and result["pyvista_open_edges"] == 0
            and edge_counts["nonmanifold_edges"] == 0
            and result["connected_components"] == 1
            and result["duplicate_face_groups"] == 0
            and intersections["intersecting_triangle_count"] == 0
            and quality["degenerate_faces_area_le_1e-16"] == 0
        )
        lines.append(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                name,
                result["trimesh_is_watertight"],
                result["pyvista_open_edges"],
                edge_counts["nonmanifold_edges"],
                result["connected_components"],
                result["duplicate_face_groups"],
                intersections["intersecting_triangle_count"],
                quality["shape_quality"]["min"],
                quality["sliver_faces_quality_lt_1e-4"],
                "基础几何 QA 通过" if passed else "需要修复/复核",
            )
        )

    lines.extend(["", "## 详细结果", ""])
    for name, result in results["cases"].items():
        quality = result["triangle_quality"]
        edge_counts = result["edge_counts"]
        intersections = result["self_intersections"]
        lines.extend(
            [
                f"### `{name}`",
                "",
                f"- STL: `{result['path']}`",
                f"- faces / vertices: `{result['faces']}` / `{result['vertices']}`",
                f"- watertight: `{result['trimesh_is_watertight']}`; winding consistent: `{result['trimesh_is_winding_consistent']}`; is_volume: `{result['is_volume']}`",
                f"- open edges: `{result['pyvista_open_edges']}`; nonmanifold edges: `{edge_counts['nonmanifold_edges']}`; max edge face count: `{edge_counts['max_edge_face_count']}`",
                f"- connected components: `{result['connected_components']}`; component faces: `{result['component_faces']}`",
                f"- duplicate face groups: `{result['duplicate_face_groups']}`",
                f"- self-intersection method: `{intersections['method']}`; intersecting triangles: `{intersections['intersecting_triangle_count']}`; elapsed: `{intersections['elapsed_sec']}s`",
                f"- triangle area min / median / max: `{quality['area']['min']}` / `{quality['area']['median']}` / `{quality['area']['max']}`",
                f"- edge length min / median / max: `{quality['edge_length']['min']}` / `{quality['edge_length']['median']}` / `{quality['edge_length']['max']}`",
                f"- aspect ratio p95 / max: `{quality['aspect_ratio']['p95']}` / `{quality['aspect_ratio']['max']}`",
                f"- shape quality min / p01 / median: `{quality['shape_quality']['min']}` / `{quality['shape_quality']['p01']}` / `{quality['shape_quality']['median']}`",
                f"- degenerate faces area <= 1e-16: `{quality['degenerate_faces_area_le_1e-16']}`",
                f"- sliver faces quality < 1e-4 / < 1e-3: `{quality['sliver_faces_quality_lt_1e-4']}` / `{quality['sliver_faces_quality_lt_1e-3']}`",
                "",
            ]
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {"cases": {}, "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    for name, path in STL_CASES.items():
        print(f"[start] {name}: {path}", flush=True)
        results["cases"][name] = qa_one(name, path)
        print(f"[done] {name}: {results['cases'][name]['elapsed_sec']}s", flush=True)

    json_path = OUT_DIR / "strict_geometry_qa.json"
    md_path = OUT_DIR / "strict_geometry_qa.md"
    json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_markdown(results, md_path)
    print(f"[done] wrote {json_path.relative_to(REPO_ROOT)}", flush=True)
    print(f"[done] wrote {md_path.relative_to(REPO_ROOT)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
