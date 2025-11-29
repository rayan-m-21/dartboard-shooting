"""
This python file provides functions needed for computing expected values on a
grid containg the dartboard. It also has a function for plotting these expected
values as a heatmap.

Example usage of this is shown in grid_plots.ipynb
"""

import expected_value
import glob
import re
import csv
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.interpolate import griddata
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def worker(answers_dir_name, triangles_dir, std_dev, points_chunk, batch_id):
    """
    The work that one CPU needs to do
    """
    filename = answers_dir_name + f"/batch_{batch_id}.csv"

    total = len(points_chunk)
    next_milestone = 0

    with open(filename, "w") as f:
        writer = csv.writer(f)
        for i, (x, y) in enumerate(points_chunk):
            val = expected_value.expected_value(triangles_dir=triangles_dir, p=(x, y), std_dev=std_dev)
            writer.writerow([x, y, val])

            progress = (i + 1) / total
            if progress >= next_milestone:
                print(f"[Worker {batch_id}] {round(next_milestone * 100)}% complete")
                f.flush()
                next_milestone += 0.1

def compute_on_grid(answers_dir_name, triangles_dir, std_dev, resolution):
    """
    Stores values in answers_dir_name/grid_data.csv
    Takes ~ 1 min for 300 points to be computed

    Resolution is the number of evaluation points along each row/col
    """
    # make directory for answers to be saved
    if os.path.isdir(answers_dir_name):
        shutil.rmtree(answers_dir_name)

    answer_path = Path(answers_dir_name)

    try:
        answer_path.mkdir(parents=True, exist_ok=True)
    except Exception:
        print(Exception.args)
        return

    x_range = y_range = np.linspace(-230, 230, resolution)
    X, Y = np.meshgrid(x_range, y_range)
    points = list(zip(X.ravel(), Y.ravel()))
    print(f'num points: {len(points)}')
    print(f'should take ~ {round((len(points) / 53361) * 15480)} seconds!')

    # separate into separate CPU's
    num_workers = os.cpu_count()
    points_chunks = [points[i::num_workers] for i in range(num_workers)]

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i, chunk in enumerate(points_chunks):
            futures.append(executor.submit(worker, answers_dir_name, triangles_dir, std_dev, chunk, i))
        # Wait for all futures to complete
        for future in futures:
            future.result()

    # Compile all batches into one file
    with open(f'{answers_dir_name}/grid_data.csv', 'w') as outf:
        for file in glob.glob(f'{answers_dir_name}/batch_*'):
            with open(file, 'r') as inf:
                for line in inf:
                    print(line.strip(),file=outf)

    for file in glob.glob(f'{answers_dir_name}/batch_*'):
        os.remove(file)

def plot_heatmap(path_to_grid_data, title=None, vmin=None, vmax=None):
    # 1. Load the CSV
    df = pd.read_csv(path_to_grid_data, header=None, names=["x", "y", "score"])

    # 2. Extract x, y, and score
    x = df["x"].values
    y = df["y"].values
    z = df["score"].values

    # 3. Create a grid to interpolate onto
    xi = np.linspace(x.min(), x.max(), 461)
    yi = np.linspace(y.min(), y.max(), 461)
    xi, yi = np.meshgrid(xi, yi)

    # 4. Interpolate score values onto grid
    zi = griddata((x, y), z, (xi, yi), method='linear')

    # 5. Plot the heatmap
    if title == None:
        title=path_to_grid_data

    plt.figure(figsize=(10, 8))
    heatmap = plt.imshow(zi, extent=[x.min(), x.max(), y.min(), y.max()],
                        origin='lower', aspect='auto', cmap='plasma',
                        vmin=vmin, vmax=vmax)

    plt.colorbar(heatmap, label="Score")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.tight_layout()
    plt.show()
