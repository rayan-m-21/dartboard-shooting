import triangulate_helper as triangulate_helper
import params as params
import re
import sys
import integrate_triangle as integrate_triangle
import numpy as np
import time
from concurrent.futures import ProcessPoolExecutor
import os
import shutil
import csv
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import subprocess

POINTS_DIR = 'BoardPoints'
TRIANGLES_DIR = 'Triangles'

def get_prob_of_hitting_region_py(triangles_dir, p, std_dev, region):
    prob_of_hitting_region = 0

    # compute the integral on each triangle of the polygon
    with open(f'{triangles_dir}/{region}.triangles', 'r') as f:
        for line in f:
            if m := re.search(r'^([^| ]*) ([^| ]*)\|([^| ]*) ([^| ]*)\|([^| ]*) ([^| ]*)$', line.strip()):
                t1=(float(m.group(1)), float(m.group(2)))
                t2=(float(m.group(3)), float(m.group(4)))
                t3=(float(m.group(5)), float(m.group(6)))

                prob_of_hitting_region += integrate_triangle.compute_integral(p, std_dev, t1, t2, t3)
            else:
                print(f'{line.strip()} NOT recognized')
                sys.exit(1)
    return float(prob_of_hitting_region)

def expected_value(triangles_dir, p=(0,0), std_dev=1):
    sum = 0

    for region in params.regions:
            sum += get_prob_of_hitting_region_py(triangles_dir, p, std_dev, region) * params.region_value[region]

    return sum
