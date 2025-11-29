import params as params
import math
import numpy as np
import re
import sys
import shutil
import display as display
from pathlib import Path
import shutil

def get_region_contour(region_key, points_per_mm):
    """
    For a specified region, this function will return a list of the contour of
    that region in anticlockwise order.
    """
    if region_key not in params.regions:
        print(f'invalid region {region_key}')
        return

    def get_points_mm_arc(start_angle, inner_rad, outer_rad):
        ret = []

        num_out_arc=math.ceil(points_per_mm*(2*np.pi*outer_rad/20))
        num_in_arc=math.ceil(points_per_mm*(2*np.pi*inner_rad/20))
        num_straight=math.ceil(points_per_mm*(outer_rad - inner_rad))

        out_arc = [] # anticlockwise
        for i in range(num_out_arc):
            angle=start_angle + i * (18 / (num_out_arc - 1))
            out_arc.append((outer_rad * np.cos(np.deg2rad(angle)), outer_rad * np.sin(np.deg2rad(angle))))

        in_arc = [] # anticlockwise
        for i in range(num_in_arc):
            angle=start_angle + i * (18 / (num_in_arc - 1))
            in_arc.append((inner_rad * np.cos(np.deg2rad(angle)), inner_rad * np.sin(np.deg2rad(angle))))

        start_ang_line = [] # going out
        end_ang_line = [] # going out
        for i in range(1, num_straight + 1):
            r = inner_rad + i * (outer_rad-inner_rad) / (num_straight + 1)
            start_ang_line.append((r * np.cos(np.deg2rad(start_angle)),r * np.sin(np.deg2rad(start_angle))))
            end_ang_line.append((r * np.cos(np.deg2rad(start_angle + 18)),r * np.sin(np.deg2rad(start_angle + 18))))

        # we must return in a sweeping order
        ret += start_ang_line
        ret += out_arc
        ret += reversed(end_ang_line)
        ret += reversed(in_arc)
        return ret

    if region_key == 'inner_bull':
        ret = []
        m = math.ceil(points_per_mm * 40)

        for i in range(m):
            angle = (360  * i) / m
            ret.append((params.radii_mm['inner_bull']*np.cos(np.deg2rad(angle)), params.radii_mm['inner_bull']*np.sin(np.deg2rad(angle))))
        return ret
    elif region_key == 'outer_bull_L':
        ret = []
        num_vert=math.ceil(points_per_mm * 9.55)
        num_in_circle=math.ceil(points_per_mm * 20)
        num_out_circle=math.ceil(points_per_mm * 50)

        out_circle = [] # anticlockwise
        for i in range(num_out_circle):
            angle=90 + i * (180 / (num_out_circle - 1))
            out_circle.append((params.radii_mm['outer_bull'] * np.cos(np.deg2rad(angle)), params.radii_mm['outer_bull'] * np.sin(np.deg2rad(angle))))

        in_circle = [] # anticlockwise
        for i in range(num_in_circle):
            angle=90 + i * (180 / (num_in_circle - 1))
            in_circle.append((params.radii_mm['inner_bull'] * np.cos(np.deg2rad(angle)), params.radii_mm['inner_bull'] * np.sin(np.deg2rad(angle))))

        top_vert = [] # down to up
        bottom_vert = [] # up to down
        for i in range(1, num_vert + 1):
            y=params.radii_mm['inner_bull'] + (params.radii_mm['outer_bull'] - params.radii_mm['inner_bull']) * i / (num_vert + 1)
            top_vert.append((0 , y))
            bottom_vert.append((0, -y))

        # we must return in a sweeping order
        ret += top_vert
        ret += out_circle
        ret += reversed(bottom_vert)
        ret += reversed(in_circle)

        return ret
    elif region_key == 'outer_bull_R':
        return list(map(lambda x: (-x[0], -x[1]), get_region_contour('outer_bull_L', points_per_mm)))
    elif m := re.search(r'double_(\d*)', region_key):
        start_angle=params.start_angles[int(m.group(1))]
        return get_points_mm_arc(start_angle, params.radii_mm['double_inner'], params.radii_mm['double_outer'])
    elif m := re.search(r'triple_(\d*)', region_key):
        start_angle=params.start_angles[int(m.group(1))]
        return get_points_mm_arc(start_angle, params.radii_mm['triple_inner'], params.radii_mm['triple_outer'])
    elif m := re.search(r'inner_(\d*)', region_key):
        start_angle=params.start_angles[int(m.group(1))]
        return get_points_mm_arc(start_angle, params.radii_mm['outer_bull'], params.radii_mm['triple_inner'])
    elif m := re.search(r'outer_(\d*)', region_key):
        start_angle=params.start_angles[int(m.group(1))]
        return get_points_mm_arc(start_angle, params.radii_mm['triple_outer'], params.radii_mm['double_inner'])
    else:
        print(f'have not coded for region key {region_key}')

def output_all_contours(points_per_mm, dir_to_write=None):
    """
    Given a value for points_per_mm to put on the diagram (must be >= 0.4),
    this function will output points along the contours of each regions
    going in an anticlockwise direction.

    By default it will save results in the directory
    output/region_contours/points_per_mm=0.8
    Each reg REG will have its points saved in a file REG.points
    - points are written individually on each lines with spaces between them
    - the points are represented with cartesian coords x,y in mm corresponding
    to offset from the center of the board

    Returns: <The path to the directory where the contours are saved>
    """
    if (points_per_mm < 0.4):
        print('Error: We need more points_per_mm for a good plot')
        sys.exit(1)

    if dir_to_write is None:
        dir_to_write = f'output/region_contours/points_per_mm={str(points_per_mm)}'

    output_path = Path(dir_to_write)

    # Remove the directory if it exists
    if output_path.exists() and output_path.is_dir():
        shutil.rmtree(output_path)

    # Create the directory (including parents if needed)
    output_path.mkdir(parents=True, exist_ok=True)

    for reg in params.regions:
        with open(str(output_path) + '/' + reg + '.points', 'w') as f:
            for (x, y) in get_region_contour(reg, points_per_mm):
                print(f'{float(x)} {float(y)}', file=f)

    return str(output_path)

if __name__ == "__main__":
    # For testing the points_per_mm visually
    i = display.get_dartboard_image()
    display.draw_sequence(i, list(map(display.mm_to_pix, get_region_contour('double_20', 0.8))))
    import cv2
    cv2.imshow('double_20', i)
    cv2.waitKey(0)
