import os
import glob
import re
import xml.etree.ElementTree as ET

import pyvista as pv


def _extract_time_from_mesh(path):
    """
    Return the time value stored in the mesh field data if available.
    """
    try:
        mesh = pv.read(path)
    except Exception:
        return None
    try:
        if "time" in mesh.field_data:
            values = mesh.field_data["time"]
            # Handle scalar or 1-element arrays
            if hasattr(values, "shape") and values.size > 0:
                return float(values.flat[0])
            return float(values)
    except Exception:
        return None
    return None


def _extract_index_from_name(path):
    """
    Extract integer index from a file name like ``time_point_12.vtp``.
    """
    name = os.path.basename(path)
    match = re.search(r"time_point_(\d+)\.vtp$", name)
    if not match:
        return 0
    return int(match.group(1))


def main():
    base_dir = os.getcwd()
    timeseries_dir = os.path.join(base_dir, "timeseries")

    if not os.path.isdir(timeseries_dir):
        raise FileNotFoundError(
            "No 'timeseries' directory found. Run 'plot_0d_results_to_3d.py' "
            "to generate time_point_*.vtp files first."
        )

    pattern = os.path.join(timeseries_dir, "time_point_*.vtp")
    files = sorted(glob.glob(pattern), key=_extract_index_from_name)

    if not files:
        raise FileNotFoundError(
            "No 'time_point_*.vtp' files found in the 'timeseries' directory. "
            "Run 'plot_0d_results_to_3d.py' first."
        )

    # Build a VTK collection file (.pvd) that references each time point
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = ET.SubElement(root, "Collection")

    for path in files:
        time_value = _extract_time_from_mesh(path)
        if time_value is None:
            # Fall back to using the index when no explicit time is stored
            time_value = _extract_index_from_name(path)
        ET.SubElement(
            collection,
            "DataSet",
            timestep=str(time_value),
            part="0",
            file=os.path.basename(path),
        )

    tree = ET.ElementTree(root)
    out_path = os.path.join(timeseries_dir, "timeseries.pvd")
    tree.write(out_path, xml_declaration=True, encoding="UTF-8")
    print(f"Created '{out_path}' with {len(files)} time steps.")


if __name__ == "__main__":
    main()
