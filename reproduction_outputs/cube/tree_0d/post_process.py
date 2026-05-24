import pyvista as pv
import numpy as np
import os
from tqdm import tqdm
import csv
from collections import defaultdict, OrderedDict


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

# Storing 5 data results: time, pressure_in, pressure_out, flow_in, flow_out
results = np.zeros((geom_data.shape[0], geom_data.shape[1] + 5, len(data['time'][0][0])))
for idx in tqdm(range(len(data['time'][0][0])),desc="Building Timeseries ",position=0):
    time = data['time'][0][0][idx]
    # Assign data from the csv file to the vessels in the results array for a given time point
    for jdx in range(results.shape[0]):
        results[jdx, 0:3, idx] = geom_data[jdx, 0:3]
        results[jdx, 3:6, idx] = geom_data[jdx, 3:6]
        results[jdx, 6, idx] = geom_data[jdx, 6]
        results[jdx, 7, idx] = geom_data[jdx, 7]
        results[jdx, 8, idx] = data['time'][jdx][0][idx]
        results[jdx, 9, idx] = data['pressure_in'][jdx][0][idx]
        results[jdx, 10, idx] = data['pressure_out'][jdx][0][idx]
        results[jdx, 11, idx] = data['flow_in'][jdx][0][idx]
        results[jdx, 12, idx] = data['flow_out'][jdx][0][idx]

# Save the results array to a numpy file
np.save("results.npy", results)
