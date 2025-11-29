import params as params
import math
import numpy as np
import re
import os
import params as params
import sys
import glob
import shutil
import subprocess
import display as display
from pathlib import Path

# Given a directory where the points of the regions are stored
# This function will generate the traiangles defined by the polygon
# and put the answers inside a specified directory
# Ant region R has its triangles saved in a file called R.triangles
def output_triangles(points_dir, triangle_points_dir):
    """
    Given a directory where all contour polygon points are stored for each region
    This function will call the binary executables to compute the triangles that
    compose each of the regions. These triangles divide each region with no
    overlap.

    Outputs are by default stored in a dir called:
    Triangles/double_2.triangles/
    Some region REG will have its triangles stored in a file called
    REG.triangles

    Returns: <path to dir where triangle points are stored>
    """
    if not os.path.isdir(points_dir):
        print(f'{points_dir} is not a directory')
        sys.exit(1)

    output_path = Path(triangle_points_dir)

    # Remove the directory if it exists
    if output_path.exists() and output_path.is_dir():
        shutil.rmtree(output_path)

    # Create the directory (including parents if needed)
    output_path.mkdir(parents=True, exist_ok=True)

    for points_filename in glob.glob(points_dir + '/*'):
        if m := re.search(rf'{points_dir}/(.*)\.points', points_filename):
            region_name = m.group(1)

            if region_name not in params.regions:
                print(f'{region_name} is not a recognized region')
                sys.exit(1)

            # We now call the c++ program which outputs the triangles
            output_filename = triangle_points_dir + '/' + region_name + '.triangles'
            subprocess.run(['./get_triangles', points_filename, output_filename])
        else:
            print(f'{points_filename} is not valid')
            sys.exit(1)
    return str(output_path)

if __name__ == "__main__":
    # example usage of output_triangles after computing contour polygons
    import approximate_regions as approximate_regions
    contours_path = approximate_regions.output_all_contours(0.9)
    print(contours_path)
    triangles_path = output_triangles(contours_path)
    print(triangles_path)
