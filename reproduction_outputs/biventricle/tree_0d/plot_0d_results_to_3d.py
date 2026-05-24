import pyvista as pv
import numpy as np
import os
import sys
from tqdm import tqdm
import csv
from collections import defaultdict


def read_csv_into_nested_dict(file_path):
    # Initialize a nested dictionary structure
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    time_steps = defaultdict(lambda: defaultdict(dict))
    time_counter = defaultdict(int)
    with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['name']
            time = float(row['time'])
            flow_in = float(row['flow_in'])
            flow_out = float(row['flow_out'])
            pressure_in = float(row['pressure_in'])
            pressure_out = float(row['pressure_out'])
            parts = name.split('_')
            branch = int(parts[0].replace('branch', ''))
            segment = int(parts[1].replace('seg', ''))
            # Check if the time has already been encountered for this branch and segment
            if time not in time_steps[branch][segment]:
                time_steps[branch][segment][time] = time_counter[(branch, segment)]
                time_counter[(branch, segment)] += 1
            time_index = time_steps[branch][segment][time]
            data['time'][branch][segment][time_index] = time
            data['flow_in'][branch][segment][time_index] = flow_in
            data['flow_out'][branch][segment][time_index] = flow_out
            data['pressure_in'][branch][segment][time_index] = pressure_in
            data['pressure_out'][branch][segment][time_index] = pressure_out
    return data


geom_data = np.genfromtxt("geom.csv",delimiter=",")
data = read_csv_into_nested_dict("output.csv")

render_pngs = os.environ.get("SVV_0D_RENDER_SCREENSHOTS", "1") == "1"
write_total = os.environ.get("SVV_0D_WRITE_TOTAL_VTP", "1") == "1"
cylinder_resolution = max(3, int(os.environ.get("SVV_0D_CYLINDER_RESOLUTION", "24")))
cylinder_capping = os.environ.get("SVV_0D_CYLINDER_CAPPING", "0") == "1"
max_units = int(os.environ.get("SVV_0D_MAX_VESSELS", "0"))
max_timesteps = int(os.environ.get("SVV_0D_MAX_TIMESTEPS", "0"))
timestep_stride = max(1, int(os.environ.get("SVV_0D_TIMESTEP_STRIDE", "1")))
disable_tqdm = os.environ.get("SVV_0D_DISABLE_TQDM", "0") == "1" or not sys.stdout.isatty()

for folder in ("timeseries", "timeseries_for_pressure_gif", "timeseries_for_flow_gif", "timeseries_for_wss_gif"):
    os.makedirs(folder, exist_ok=True)

all_units = []
for vessel in sorted(data['flow_in'].keys()):
    start = geom_data[vessel, 0:3]
    end = geom_data[vessel, 3:6]
    direction = geom_data[vessel, 3:6] - geom_data[vessel, 0:3]
    direction_norm = np.linalg.norm(direction)
    if direction_norm == 0.0:
        continue
    direction = direction / direction_norm
    length = float(geom_data[vessel, 6])
    radius = float(geom_data[vessel, 7])
    center = 0.5 * direction * length + start
    for segment in sorted(data['flow_in'][vessel].keys()):
        base_mesh = pv.Cylinder(
            center=center,
            direction=direction,
            height=length,
            radius=radius,
            resolution=cylinder_resolution,
            capping=cylinder_capping,
        )
        base_mesh = base_mesh.elevation(low_point=end, high_point=start, scalar_range=[0.0, 1.0])
        axis_unit = np.asarray(base_mesh.point_data['Elevation'])
        all_units.append({
            'vessel': vessel,
            'segment': segment,
            'radius': radius,
            'base_mesh': base_mesh,
            'axis_unit': axis_unit,
        })
        if max_units > 0 and len(all_units) >= max_units:
            break
    if max_units > 0 and len(all_units) >= max_units:
        break

if len(all_units) == 0:
    raise RuntimeError("No vessel geometry could be generated from geom.csv.")

first_vessel = next(iter(data['time']))
first_segment = next(iter(data['time'][first_vessel]))
num_timesteps = len(data['time'][first_vessel][first_segment])
if max_timesteps > 0:
    num_timesteps = min(num_timesteps, max_timesteps)
timesteps = list(range(0, num_timesteps, timestep_stride))

print("Found {} vessel segments across {} timesteps".format(len(all_units), len(timesteps)))
print(
    "Options: render_pngs={}, write_total={}, cylinder_resolution={}, cylinder_capping={}, max_vessels={}, max_timesteps={}, timestep_stride={}".format(
        render_pngs,
        write_total,
        cylinder_resolution,
        cylinder_capping,
        max_units,
        max_timesteps,
        timestep_stride,
    )
)

min_pressure = np.inf
max_pressure = -np.inf
min_flow = np.inf
max_flow = -np.inf
min_wss = np.inf
max_wss = -np.inf
for idx in timesteps:
    for unit in all_units:
        vessel = unit['vessel']
        segment = unit['segment']
        radius = unit['radius']
        p_in = float(data['pressure_in'][vessel][segment][idx]) / 1333.33
        p_out = float(data['pressure_out'][vessel][segment][idx]) / 1333.33
        q_in = float(data['flow_in'][vessel][segment][idx])
        vel = (q_in / (np.pi * radius**2)) / 2.0
        re = (1.06 * 2.0 * radius * vel) / 0.04
        fd = 0.0 if re == 0.0 else 64.0 / re
        wss = vel * fd * 1.06
        min_pressure = min(min_pressure, p_in, p_out)
        max_pressure = max(max_pressure, p_in, p_out)
        min_flow = min(min_flow, q_in)
        max_flow = max(max_flow, q_in)
        min_wss = min(min_wss, wss)
        max_wss = max(max_wss, wss)

total = None
for out_idx, idx in enumerate(tqdm(timesteps, desc="Building/Saving Timeseries", position=0, disable=disable_tqdm)):
    time_value = float(data['time'][first_vessel][first_segment][idx])
    meshes = []
    for unit in all_units:
        vessel = unit['vessel']
        segment = unit['segment']
        radius = unit['radius']
        axis_unit = unit['axis_unit']
        mesh = unit['base_mesh'].copy(deep=True)
        p_in = float(data['pressure_in'][vessel][segment][idx]) / 1333.33
        p_out = float(data['pressure_out'][vessel][segment][idx]) / 1333.33
        q_in = float(data['flow_in'][vessel][segment][idx])
        mesh.point_data['Pressure [mmHg]'] = p_out + axis_unit * (p_in - p_out)
        mesh.cell_data['Flow [mL/s]'] = q_in
        vel = (q_in / (np.pi * radius**2)) / 2.0
        re = (1.06 * 2.0 * radius * vel) / 0.04
        fd = 0.0 if re == 0.0 else 64.0 / re
        mesh.cell_data['WSS [dyne/cm^2]'] = vel * fd * 1.06
        meshes.append(mesh)

    time_mesh = meshes[0].merge(meshes[1:]) if len(meshes) > 1 else meshes[0]
    time_mesh.field_data['time'] = np.array([time_value], dtype=float)
    time_mesh.save(os.path.join(os.getcwd(), "timeseries", "time_point_{}.vtp".format(out_idx)))

    if render_pngs:
        p = pv.Plotter(off_screen=True)
        p.add_mesh(time_mesh, scalars='Pressure [mmHg]', clim=[round(min_pressure, 4), round(max_pressure, 4)], cmap="coolwarm")
        p.show(auto_close=True, screenshot=os.path.join(os.getcwd(), "timeseries_for_pressure_gif", "time_point_{}.png".format(out_idx)))
        p = pv.Plotter(off_screen=True)
        p.add_mesh(time_mesh, scalars='Flow [mL/s]', clim=[round(min_flow, 4), round(max_flow, 4)], cmap="GnBu")
        p.show(auto_close=True, screenshot=os.path.join(os.getcwd(), "timeseries_for_flow_gif", "time_point_{}.png".format(out_idx)))
        p = pv.Plotter(off_screen=True)
        p.add_mesh(time_mesh, scalars='WSS [dyne/cm^2]', clim=[round(min_wss, 2), round(max_wss, 2)], cmap="coolwarm")
        p.show(auto_close=True, screenshot=os.path.join(os.getcwd(), "timeseries_for_wss_gif", "time_point_{}.png".format(out_idx)))

    if write_total:
        total = time_mesh if total is None else total.merge(time_mesh)

if write_total and total is not None:
    total.save(os.path.join(os.getcwd(), "timeseries", "total.vtp"))

