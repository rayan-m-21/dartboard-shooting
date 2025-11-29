import cv2
import numpy as np
import math
import time
import re
import os


# Regulation dartboard radii (in mm)
radii_mm = {
    'inner_bull': 6.35,
    'outer_bull': 15.9,
    'triple_inner': 107,
    'triple_outer': 115,
    'double_inner': 170,
    'double_outer': 178,
    'board_edge': 225.5
}

# Initialise regions set
regions = ['inner_bull', 'outer_bull_L', 'outer_bull_R']
for i in [6, 13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10]:
    regions.append(f'double_{i}')
    regions.append(f'triple_{i}')
    regions.append(f'inner_{i}')
    regions.append(f'outer_{i}')

# initialize starting angle for each score
start_angles = {}
for i, x in enumerate([6, 13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10]):
    start_angles[x] = -9 + 18*i

# initialise values for each region
region_value = { 'inner_bull': 50, 'outer_bull_L': 25, 'outer_bull_R': 25, 'outside': 0 }
for i in [6, 13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10]:
    region_value[f'inner_{i}'] = i
    region_value[f'outer_{i}'] = i
    region_value[f'double_{i}'] = 2 * i
    region_value[f'triple_{i}'] = 3 * i
