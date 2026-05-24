"""Execute the software-focused reproduction plan for svVascularize.

The script intentionally uses small generation sizes. It verifies the repository
pipeline end-to-end before attempting paper-scale benchmarks.
"""

from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh

from svv.domain.domain import Domain
from svv.forest.forest import Forest
from svv.simulation.simulation import Simulation
from svv.tree.tree import Tree
from svv.utils.remeshing.mmg import get_mmg_exe
from svv.utils.solvers.solver_0d import get_solver_0d_exe


REPO_ROOT = Path(__file__).resolve().parent
OUT_ROOT = REPO_ROOT / "reproduction_outputs"
REUSE_OUTPUTS = os.environ.get("SVV_REPRO_REUSE_OUTPUTS") == "1"


def log(message: str) -> None:
    print(message, flush=True)


def write_summary(summary: dict) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary_path = OUT_ROOT / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def timed(summary: dict, name: str):
    class _Timer:
        def __enter__(self):
            self.start = time.perf_counter()
            log(f"[start] {name}")
            return self

        def __exit__(self, exc_type, exc, _tb):
            elapsed = time.perf_counter() - self.start
            summary.setdefault("timings_sec", {})[name] = round(elapsed, 3)
            if exc is None:
                log(f"[done] {name}: {elapsed:.2f}s")
            else:
                summary.setdefault("errors", {})[name] = repr(exc)
                log(f"[fail] {name}: {exc!r}")
            return False

    return _Timer()


def make_domain(mesh: pv.PolyData, seed: int, resolution: int, summary: dict, name: str) -> Domain:
    domain = Domain(mesh.triangulate())
    domain.set_random_seed(seed)
    with timed(summary, f"{name}.domain.create"):
        domain.create()
    with timed(summary, f"{name}.domain.solve"):
        domain.solve()
    with timed(summary, f"{name}.domain.build"):
        domain.build(resolution=resolution)
    summary[name]["domain"] = {
        "points": int(domain.points.shape[0]),
        "patches": int(len(domain.patches)),
        "boundary_cells": int(domain.boundary.n_cells),
        "mesh_cells": int(domain.mesh.n_cells),
        "convexity": None if domain.convexity is None else float(domain.convexity),
    }
    return domain


def make_tree(domain: Domain, adds: int, seed: int, summary: dict, name: str) -> Tree:
    random.seed(seed)
    np.random.seed(seed)
    tree = Tree(preallocation_step=max(128, 2 * adds + 8))
    tree.random_seed = seed
    tree.set_domain(domain)
    tree.set_root()
    with timed(summary, f"{name}.tree.n_add_{adds}"):
        tree.n_add(adds)
    summary[name]["tree"] = {
        "adds": int(adds),
        "segments": int(tree.data.shape[0]),
        "terminals": int(tree.n_terminals),
        "root_flow": float(tree.data[0, 22]),
        "root_radius": float(tree.data[0, 21]),
    }
    return tree


def export_centerlines(tree: Tree, outdir: Path, summary: dict, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    centerlines, _ = tree.export_centerlines(points_per_unit_length=10)
    path = outdir / "centerlines.vtp"
    centerlines.save(path)
    summary[name]["centerlines"] = {
        "path": str(path.relative_to(REPO_ROOT)),
        "points": int(centerlines.n_points),
        "cells": int(centerlines.n_cells),
    }


def export_and_run_0d(tree: Tree, outdir: Path, summary: dict, name: str) -> None:
    case_dir = outdir / "tree_0d"
    if REUSE_OUTPUTS and (case_dir / "output.csv").is_file():
        summary[name]["zerod"] = {
            "case_dir": str(case_dir.relative_to(REPO_ROOT)),
            "returncode": 0,
            "solver_input": str((case_dir / "solver_0d.in").relative_to(REPO_ROOT)),
            "geom": str((case_dir / "geom.csv").relative_to(REPO_ROOT)),
            "output_csv": str((case_dir / "output.csv").relative_to(REPO_ROOT)),
            "output_csv_exists": True,
            "reused": True,
        }
        log(f"[skip] {name}.0d reused existing output.csv")
        return

    sim = Simulation(tree, name="tree_0d", directory=str(outdir))
    with timed(summary, f"{name}.0d.export"):
        case_dir = Path(
            sim.export_0d_fluid_simulation(
                get_0d_solver=True,
                number_cardiac_cycles=1,
                number_time_pts_per_cycle=5,
            )
        )
    with timed(summary, f"{name}.0d.run"):
        run = subprocess.run(
            [sys.executable, "run.py"],
            cwd=case_dir,
            text=True,
            capture_output=True,
            timeout=120,
            check=False,
        )
    summary[name]["zerod"] = {
        "case_dir": str(case_dir.relative_to(REPO_ROOT)),
        "returncode": int(run.returncode),
        "solver_input": str((case_dir / "solver_0d.in").relative_to(REPO_ROOT)),
        "geom": str((case_dir / "geom.csv").relative_to(REPO_ROOT)),
        "output_csv": str((case_dir / "output.csv").relative_to(REPO_ROOT)),
        "output_csv_exists": (case_dir / "output.csv").is_file(),
        "stdout_tail": run.stdout[-500:],
        "stderr_tail": run.stderr[-500:],
    }
    if run.returncode != 0:
        raise RuntimeError(f"0D solver failed for {name}: returncode={run.returncode}")


def _flatten_meshes(meshes):
    flat = []
    for item in meshes:
        if isinstance(item, (list, tuple)):
            flat.extend(_flatten_meshes(item))
        elif item is not None:
            flat.append(item)
    return flat


def export_3d(tree: Tree, outdir: Path, summary: dict, name: str) -> None:
    case_dir = outdir / "tree_3d"
    mesh_dir = case_dir / "mesh"
    mesh_path = mesh_dir / "fluid_msh_0" / "fluid_msh_0.vtu"
    xml_files = sorted(case_dir.glob("fluid_simulation_*.xml"))
    if REUSE_OUTPUTS and xml_files and mesh_path.is_file():
        mesh = pv.read(mesh_path)
        try:
            surface_mesh_cells = int(mesh.extract_surface().n_cells)
        except Exception:
            surface_mesh_cells = None
        summary[name]["threed"] = {
            "xml_files": [str(path.relative_to(REPO_ROOT)) for path in xml_files],
            "mesh_dir": str(mesh_dir.relative_to(REPO_ROOT)),
            "volume_mesh_cells": int(mesh.n_cells),
            "surface_mesh_cells": surface_mesh_cells,
            "reused": True,
        }
        log(f"[skip] {name}.3d reused existing mesh/XML")
        return

    sim = Simulation(tree, name="tree_3d", directory=str(outdir))
    with timed(summary, f"{name}.3d.export"):
        files = sim.export_3d_fluid_simulation(
            boundary_layer=False,
            remesh_vol=False,
            minratio=1.1,
            mindihedral=10.0,
        )
    mesh_dir = outdir / "tree_3d" / "mesh"
    summary[name]["threed"] = {
        "xml_files": [str(Path(p).relative_to(REPO_ROOT)) for p in files],
        "mesh_dir": str(mesh_dir.relative_to(REPO_ROOT)),
        "volume_mesh_cells": int(sim.fluid_domain_volume_meshes[0].n_cells),
        "surface_mesh_cells": int(sim.fluid_domain_surface_meshes[0].n_cells),
    }


def export_stl_candidate(tree: Tree, outdir: Path, summary: dict, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    stl_path = outdir / "candidate_vasculature.stl"
    if REUSE_OUTPUTS and stl_path.is_file():
        surface = pv.read(stl_path)
        log(f"[skip] {name}.stl_candidate reused existing STL")
    else:
        with timed(summary, f"{name}.stl_candidate.export"):
            surface = tree.export_solid(watertight=True)
        surface.save(stl_path)

    tri = trimesh.load_mesh(str(stl_path), force="mesh")
    connected = tri.split(only_watertight=False)
    try:
        trimesh.repair.fix_normals(tri)
        self_intersection_check = "not_available"
        has_self_intersections = None
    except Exception as exc:
        self_intersection_check = f"failed: {exc!r}"
        has_self_intersections = None

    radii = np.asarray(tree.data[:, 21], dtype=float)
    summary[name]["stl_candidate"] = {
        "path": str(stl_path.relative_to(REPO_ROOT)),
        "surface_cells": int(surface.n_cells),
        "surface_points": int(surface.n_points),
        "open_edges": int(surface.n_open_edges),
        "pyvista_is_watertight": int(surface.n_open_edges) == 0,
        "trimesh_is_watertight": bool(tri.is_watertight),
        "connected_components": int(len(connected)),
        "min_radius": float(np.min(radii)),
        "min_diameter": float(2.0 * np.min(radii)),
        "self_intersection_check": self_intersection_check,
        "has_self_intersections": has_self_intersections,
    }


def summarize_stl_candidate(stl_path: Path, surface: pv.PolyData, radii: np.ndarray) -> dict:
    tri = trimesh.load_mesh(str(stl_path), force="mesh")
    connected = tri.split(only_watertight=False)
    try:
        trimesh.repair.fix_normals(tri)
        self_intersection_check = "not_available"
        has_self_intersections = None
    except Exception as exc:
        self_intersection_check = f"failed: {exc!r}"
        has_self_intersections = None

    radii = np.asarray(radii, dtype=float)
    radii = radii[np.isfinite(radii)]
    return {
        "path": str(stl_path.relative_to(REPO_ROOT)),
        "surface_cells": int(surface.n_cells),
        "surface_points": int(surface.n_points),
        "open_edges": int(surface.n_open_edges),
        "pyvista_is_watertight": int(surface.n_open_edges) == 0,
        "trimesh_is_watertight": bool(tri.is_watertight),
        "connected_components": int(len(connected)),
        "min_radius": None if radii.size == 0 else float(np.min(radii)),
        "min_diameter": None if radii.size == 0 else float(2.0 * np.min(radii)),
        "self_intersection_check": self_intersection_check,
        "has_self_intersections": has_self_intersections,
    }


def export_forest_stl_candidate(forest: Forest, outdir: Path, summary: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    stl_path = outdir / "candidate_vasculature.stl"
    if REUSE_OUTPUTS and stl_path.is_file():
        surface = pv.read(stl_path)
        log("[skip] cube_forest.stl_candidate reused existing STL")
    else:
        with timed(summary, "cube_forest.stl_candidate.export"):
            solids, _, _ = forest.connections.export_solid(extrude_roots=False)
            solids = [solid for solid in solids if solid is not None]
            if not solids:
                raise RuntimeError("No forest solids were generated for STL export.")
            surface = solids[0] if len(solids) == 1 else pv.merge(solids)
            surface = surface.triangulate().compute_normals(auto_orient_normals=True)
        surface.save(stl_path)

    radii = []
    for network in forest.networks:
        for tree in network:
            radii.extend(np.asarray(tree.data[:, 21], dtype=float).tolist())
    summary["cube_forest"]["stl_candidate"] = summarize_stl_candidate(stl_path, surface, np.asarray(radii))


def make_forest(domain: Domain, summary: dict) -> Forest:
    forest = Forest(
        domain=domain,
        n_networks=1,
        n_trees_per_network=[2],
        preallocation_step=128,
    )
    forest.set_domain(domain)
    with timed(summary, "cube_forest.roots"):
        forest.set_roots()
    with timed(summary, "cube_forest.add_1"):
        forest.add(1)
    with timed(summary, "cube_forest.connect"):
        forest.connect()
    summary["cube_forest"] = {
        "networks": int(forest.n_networks),
        "trees_per_network": [int(v) for v in forest.n_trees_per_network],
        "segments_per_tree": [int(tree.data.shape[0]) for tree in forest.networks[0]],
        "has_connections": forest.connections is not None,
        "connection_count": int(len(forest.connections.tree_connections)),
    }
    return forest


def export_forest_0d(forest: Forest, outdir: Path, summary: dict) -> None:
    sim = Simulation(forest, name="forest_0d", directory=str(outdir))
    with timed(summary, "cube_forest.0d.export"):
        case_dir = Path(
            sim.export_0d_fluid_simulation(
                0,
                [0],
                get_0d_solver=True,
                number_cardiac_cycles=1,
                number_time_pts_per_cycle=5,
            )
        )
    summary["cube_forest"]["zerod"] = {
        "case_dir": str(case_dir.relative_to(REPO_ROOT)),
        "solver_input": str((case_dir / "solver_0d.in").relative_to(REPO_ROOT)),
        "geom": str((case_dir / "geom.csv").relative_to(REPO_ROOT)),
        "output_csv": str((case_dir / "output.csv").relative_to(REPO_ROOT)),
        "solver_input_exists": (case_dir / "solver_0d.in").is_file(),
        "output_csv_exists": (case_dir / "output.csv").is_file(),
    }
    with timed(summary, "cube_forest.0d.run"):
        run = subprocess.run(
            [sys.executable, "run.py"],
            cwd=case_dir,
            text=True,
            capture_output=True,
            timeout=120,
            check=False,
        )
    summary["cube_forest"]["zerod"].update(
        {
            "status": "success" if run.returncode == 0 else "failed",
            "returncode": int(run.returncode),
            "output_csv_exists": (case_dir / "output.csv").is_file(),
            "stdout_tail": run.stdout[-500:],
            "stderr_tail": run.stderr[-500:],
        }
    )
    if run.returncode != 0:
        log(f"[warn] cube_forest.0d.run failed: returncode={run.returncode}")


def export_forest_3d_attempt(forest: Forest, outdir: Path, summary: dict) -> None:
    sim = Simulation(forest, name="forest_3d", directory=str(outdir))
    mesh_dir = outdir / "forest_3d" / "mesh"
    mesh_files = sorted(mesh_dir.glob("fluid_msh_*/*.vtu"))
    xml_files = sorted((outdir / "forest_3d").glob("fluid_simulation_*.xml"))
    if REUSE_OUTPUTS and xml_files and mesh_files:
        volume_meshes = [pv.read(path) for path in mesh_files]
        surface_mesh_cells = []
        for mesh in volume_meshes:
            try:
                surface_mesh_cells.append(int(mesh.extract_surface().n_cells))
            except Exception:
                surface_mesh_cells.append(None)
        summary["cube_forest"]["threed_attempt"] = {
            "status": "success",
            "xml_files": [str(path.relative_to(REPO_ROOT)) for path in xml_files],
            "mesh_dir": str(mesh_dir.relative_to(REPO_ROOT)),
            "volume_mesh_cells": [int(mesh.n_cells) for mesh in volume_meshes],
            "surface_mesh_cells": surface_mesh_cells,
            "reused": True,
        }
        log("[skip] cube_forest.3d reused existing mesh/XML")
        return
    try:
        with timed(summary, "cube_forest.3d.export_attempt"):
            files = sim.export_3d_fluid_simulation(
                boundary_layer=False,
                remesh_vol=False,
                minratio=1.1,
                mindihedral=10.0,
            )
        volume_meshes = _flatten_meshes(sim.fluid_domain_volume_meshes)
        surface_meshes = _flatten_meshes(sim.fluid_domain_surface_meshes)
        summary["cube_forest"]["threed_attempt"] = {
            "status": "success",
            "xml_files": [str(Path(p).relative_to(REPO_ROOT)) for p in files],
            "mesh_dir": str((outdir / "forest_3d" / "mesh").relative_to(REPO_ROOT)),
            "volume_mesh_cells": [int(mesh.n_cells) for mesh in volume_meshes],
            "surface_mesh_cells": [int(mesh.n_cells) for mesh in surface_meshes],
        }
    except Exception as exc:
        summary.setdefault("errors", {})["cube_forest.3d.export_attempt"] = repr(exc)
        summary["cube_forest"]["threed_attempt"] = {
            "status": "failed",
            "error": repr(exc),
            "note": "Connected forest 3D export is treated as best effort in this reproduction stage.",
        }
        log(f"[warn] cube_forest.3d.export_attempt failed: {exc!r}")


def main() -> int:
    os.environ.setdefault("SVV_TELEMETRY_DISABLED", "1")
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    if OUT_ROOT.exists() and not REUSE_OUTPUTS:
        shutil.rmtree(OUT_ROOT)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    summary: dict = {
        "environment": {
            "python": sys.version,
            "pyvista": pv.__version__,
            "mmg": {tool: str(get_mmg_exe(tool).relative_to(REPO_ROOT)) for tool in ("mmg2d", "mmg3d", "mmgs")},
            "solver_0d": str(get_solver_0d_exe().relative_to(REPO_ROOT)),
        },
        "cube": {},
        "cube_3d": {},
        "cube_forest": {},
        "biventricle": {},
        "stl_export_check": {},
    }

    cube_domain = make_domain(
        pv.Cube().triangulate(),
        seed=101,
        resolution=15,
        summary=summary,
        name="cube",
    )
    cube_tree = make_tree(cube_domain, adds=3, seed=102, summary=summary, name="cube")
    export_centerlines(cube_tree, OUT_ROOT / "cube" / "centerlines", summary, "cube")
    export_and_run_0d(cube_tree, OUT_ROOT / "cube", summary, "cube")
    export_3d(cube_tree, OUT_ROOT / "cube", summary, "cube")
    export_stl_candidate(cube_tree, OUT_ROOT / "stl_export_check" / "cube", summary, "stl_export_check")
    write_summary(summary)

    cube_3d_tree = make_tree(cube_domain, adds=1, seed=103, summary=summary, name="cube_3d")
    export_centerlines(cube_3d_tree, OUT_ROOT / "cube_3d" / "centerlines", summary, "cube_3d")
    export_and_run_0d(cube_3d_tree, OUT_ROOT / "cube_3d", summary, "cube_3d")
    export_3d(cube_3d_tree, OUT_ROOT / "cube_3d", summary, "cube_3d")
    write_summary(summary)

    forest = make_forest(cube_domain, summary)
    export_forest_0d(forest, OUT_ROOT / "cube_forest", summary)
    export_forest_3d_attempt(forest, OUT_ROOT / "cube_forest", summary)
    export_forest_stl_candidate(forest, OUT_ROOT / "stl_export_check" / "cube_forest", summary)
    write_summary(summary)

    biv_mesh = pv.read(str(REPO_ROOT / "biventricle_tissue.stl")).triangulate()
    biv_domain = make_domain(
        biv_mesh,
        seed=201,
        resolution=12,
        summary=summary,
        name="biventricle",
    )
    biv_tree = make_tree(biv_domain, adds=1, seed=202, summary=summary, name="biventricle")
    export_centerlines(biv_tree, OUT_ROOT / "biventricle" / "centerlines", summary, "biventricle")
    export_and_run_0d(biv_tree, OUT_ROOT / "biventricle", summary, "biventricle")
    export_3d(biv_tree, OUT_ROOT / "biventricle", summary, "biventricle")
    export_stl_candidate(biv_tree, OUT_ROOT / "stl_export_check" / "biventricle", summary, "biventricle")
    write_summary(summary)

    log(f"[done] wrote {(OUT_ROOT / 'summary.json').relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
