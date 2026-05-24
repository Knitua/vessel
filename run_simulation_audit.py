#!/usr/bin/env python3
"""Audit generated 0D results and 3D svMultiPhysics inputs."""

from __future__ import annotations

import csv
import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

from svv.utils.solvers.solver_0d import get_solver_0d_candidates


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "reproduction_outputs"
REPORT_DIR = OUT / "simulation_audit"


CASES_0D = {
    "cube": OUT / "cube" / "tree_0d" / "output.csv",
    "cube_3d": OUT / "cube_3d" / "tree_0d" / "output.csv",
    "cube_forest": OUT / "cube_forest" / "forest_0d" / "output.csv",
    "biventricle": OUT / "biventricle" / "tree_0d" / "output.csv",
}


CASES_3D = {
    "cube": OUT / "cube" / "tree_3d" / "fluid_simulation_0-0.xml",
    "cube_3d": OUT / "cube_3d" / "tree_3d" / "fluid_simulation_0-0.xml",
    "cube_forest": OUT / "cube_forest" / "forest_3d" / "fluid_simulation_0-0.xml",
    "biventricle": OUT / "biventricle" / "tree_3d" / "fluid_simulation_0-0.xml",
}


def _range(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "max": None}
    return {"min": min(values), "max": max(values)}


def audit_0d(path: Path) -> dict:
    if not path.exists():
        return {"exists": False}
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    flows_in = [float(row["flow_in"]) for row in rows]
    flows_out = [float(row["flow_out"]) for row in rows]
    pressures_in = [float(row["pressure_in"]) for row in rows]
    pressures_out = [float(row["pressure_out"]) for row in rows]
    max_abs_balance_error = max(
        (abs(a - b) for a, b in zip(flows_in, flows_out)),
        default=0.0,
    )
    flow_scale = max((abs(v) for v in flows_in), default=0.0) or 1.0
    names = sorted({row["name"] for row in rows})
    times = sorted({float(row["time"]) for row in rows})
    return {
        "exists": True,
        "path": str(path.relative_to(ROOT)),
        "rows": len(rows),
        "vessels": len(names),
        "timepoints": len(times),
        "times": times,
        "flow_in": _range(flows_in),
        "flow_out": _range(flows_out),
        "pressure_in": _range(pressures_in),
        "pressure_out": _range(pressures_out),
        "max_abs_flow_balance_error": max_abs_balance_error,
        "relative_flow_balance_error": max_abs_balance_error / flow_scale,
    }


def audit_3d(path: Path) -> dict:
    if not path.exists():
        return {"exists": False}
    root = ET.parse(path).getroot()
    general = root.find("GeneralSimulationParameters")
    mesh = root.find("Add_mesh")
    equation = root.find("Add_equation")
    faces = [node.attrib.get("name") for node in mesh.findall("Add_face")] if mesh is not None else []
    bcs = []
    if equation is not None:
        for node in equation.findall("Add_BC"):
            bcs.append(
                {
                    "name": node.attrib.get("name"),
                    "type": node.findtext("Type"),
                    "value": node.findtext("Value"),
                }
            )
    bc_names = {bc["name"] for bc in bcs}
    mesh_rel = mesh.findtext("Mesh_file_path") if mesh is not None else None
    mesh_path = path.parent / mesh_rel if mesh_rel else None
    result_files = sorted(path.parent.glob("result*"))
    return {
        "exists": True,
        "path": str(path.relative_to(ROOT)),
        "number_of_time_steps": general.findtext("Number_of_time_steps") if general is not None else None,
        "time_step_size": general.findtext("Time_step_size") if general is not None else None,
        "mesh_name": mesh.attrib.get("name") if mesh is not None else None,
        "mesh_file": str(mesh_path.relative_to(ROOT)) if mesh_path is not None else None,
        "mesh_file_exists": mesh_path.exists() if mesh_path is not None else False,
        "faces": faces,
        "boundary_conditions": bcs,
        "missing_bc_faces": [face for face in faces if face not in bc_names],
        "has_result_files": bool(result_files),
        "result_files": [str(p.relative_to(ROOT)) for p in result_files],
    }


def audit_solvers() -> dict:
    solver_0d = get_solver_0d_candidates()
    return {
        "zero_d_candidates": [
            {"path": str(path), "exists": path.is_file()} for path in solver_0d.candidates
        ],
        "svmultiphysics_path": shutil.which("svMultiPhysics") or shutil.which("svmultiphysics"),
    }


def write_markdown(data: dict, path: Path) -> None:
    lines = [
        "# 仿真部分审计报告",
        "",
        "本报告检查 reproduction 输出中的 0D 求解结果和 3D svMultiPhysics 输入文件。",
        "",
        "结论：0D 已经实际运行并生成 `output.csv`；3D 目前已生成网格与 XML 输入文件，但本机未发现 `svMultiPhysics` 可执行程序，也未发现 3D `result*` 输出，因此当前不能声明 3D CFD 已完成。",
        "",
        "## 求解器状态",
        "",
    ]
    solver = data["solvers"]
    lines.append(f"- `svMultiPhysics`: `{solver['svmultiphysics_path'] or 'not found in PATH'}`")
    lines.append("- 0D solver candidates:")
    for candidate in solver["zero_d_candidates"]:
        lines.append(f"  - `{candidate['path']}` exists={candidate['exists']}")
    lines.extend(["", "## 0D 求解结果", ""])
    lines.append("| 体系 | rows | vessels | timepoints | flow balance rel. error | pressure_in range | pressure_out range |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for name, result in data["zero_d"].items():
        if not result.get("exists"):
            lines.append(f"| `{name}` | missing |  |  |  |  |  |")
            continue
        pin = result["pressure_in"]
        pout = result["pressure_out"]
        lines.append(
            f"| `{name}` | {result['rows']} | {result['vessels']} | {result['timepoints']} | "
            f"{result['relative_flow_balance_error']:.3e} | "
            f"{pin['min']:.6g} to {pin['max']:.6g} | {pout['min']:.6g} to {pout['max']:.6g} |"
        )
    lines.extend(["", "## 3D 输入文件", ""])
    lines.append("| 体系 | mesh exists | steps | dt | faces | missing BC faces | result files |")
    lines.append("|---|---:|---:|---:|---:|---|---:|")
    for name, result in data["three_d"].items():
        if not result.get("exists"):
            lines.append(f"| `{name}` | missing |  |  |  |  |  |")
            continue
        missing = ", ".join(f"`{face}`" for face in result["missing_bc_faces"]) or "none"
        lines.append(
            f"| `{name}` | {result['mesh_file_exists']} | {result['number_of_time_steps']} | "
            f"{result['time_step_size']} | {len(result['faces'])} | {missing} | {len(result['result_files'])} |"
        )
    lines.extend(
        [
            "",
            "## 判断",
            "",
            "- 0D：已经完成一次 steady/准稳态 lumped-parameter 求解；这些结果可以用于检查网络阻力、压力和流量分配是否数量级合理。",
            "- 3D：当前处于 simulation input/export 阶段；XML、体网格和面文件已经存在，但没有本机求解器和结果文件。",
            "- 当前 XML 已补齐 `lumen_*` 血管壁面的 no-slip Dirichlet 边界条件；后续重新导出也会由代码自动补齐。",
            "- 下一步若要做真正 3D CFD，需要安装或提供 `svMultiPhysics`，然后从最小的 `cube_3d` 开始试跑 1-10 step，检查收敛、质量守恒、压力/速度场，再扩大到 `cube_forest` 和 `biventricle`。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    data = {
        "solvers": audit_solvers(),
        "zero_d": {name: audit_0d(path) for name, path in CASES_0D.items()},
        "three_d": {name: audit_3d(path) for name, path in CASES_3D.items()},
    }
    json_path = REPORT_DIR / "simulation_audit.json"
    md_path = REPORT_DIR / "simulation_audit.md"
    json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    write_markdown(data, md_path)
    print(f"Wrote {md_path.relative_to(ROOT)}")
    print(f"Wrote {json_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
