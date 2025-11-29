import params
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import os
import tempfile
import expected_value


# Function for determining what region we are in based on the (x,y) in mm
# this works by scanning through each of the .triangles files to see which triangle it belongs to

def is_in_triangle(a, b, c, point):
    def cross_prod(u, v):
        return u[0]*v[1] - u[1]*v[0]

    def get_triangle_area(p, q, t):
        u = [q[0] - p[0], q[1] - p[1]]
        v = [t[0] - p[0], t[1] - p[1]]
        return 0.5 * abs(cross_prod(u, v))

    A1 = get_triangle_area(a, b, point)
    A2 = get_triangle_area(a, c, point)
    A3 = get_triangle_area(c, b, point)
    totalA = get_triangle_area(a, b, c)

    return abs(totalA - A1 - A2 - A3) < 10**(-2)

def get_region(point, triangles_dir):
    if point[0]**2 + point[1]**2 > params.radii_mm['double_outer']**2:
        return 'outside'

    for f in glob.glob(f'{triangles_dir}/*.triangles'):
        if m := re.search(fr'{triangles_dir}/(.*)\.triangles', f):
            region = m.group(1)

            with open(f, 'r') as stream:
                for line in stream.readlines():
                    if m2 := re.search(r'^([^| ]*) ([^| ]*)\|([^| ]*) ([^| ]*)\|([^| ]*) ([^| ]*)$', line.strip()):
                        t1=(float(m2.group(1)), float(m2.group(2)))
                        t2=(float(m2.group(3)), float(m2.group(4)))
                        t3=(float(m2.group(5)), float(m2.group(6)))

                        if is_in_triangle(t1, t2, t3, point):
                            return region
                    else:
                        print('failed to read triangles file')
                        return None
        else:
            print('failed to read directory correctly')
            return None

    # This case can occur when the point is in the board but not sitting in any
    # of the triangles
    # if this is the case we should trigger a regeneration
    return None

def generate_throw(p, sigma):
    # The covariance matrix for N_2(p, sigma*I) is a diagonal matrix with
    # variance (sigma^2) on the diagonal.
    covariance_matrix = [[sigma**2, 0], [0, sigma**2]]

    # Generate a single sample (throw) from the multivariate normal distribution.
    return np.random.multivariate_normal(mean=p, cov=covariance_matrix)

def get_random_score(p, sigma):
    reg = None

    while reg is None:
        point = generate_throw(p, sigma)
        reg = get_region(point, 'output/triangles/points_per_mm=5')

    return params.region_value[reg]

def worker(args):
    p, sigma, N, fileToWrite = args

    with open(fileToWrite, "w") as f:
        for _ in range(N):
            print(get_random_score(np.array(p), sigma), file=f)


def get_N_random_scores(p, sigma, N):
    n_workers = os.cpu_count()
    temp_files = []

    # Prepare filenames only
    for _ in range(n_workers):
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_files.append(temp_file.name)  # save only the name
        temp_file.close()  # close it, worker will open it

    # Prepare arguments for each worker
    args_list = [(p, sigma, int(N / n_workers), fname) for fname in temp_files]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(worker, args) for args in args_list]
        for future in futures:
            future.result()

    scores = []
    for filename in temp_files:
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                scores.append(int(line))

    # fill in extra scores that were a result of rounding error
    for _ in range(N - len(scores)):
        scores.append(get_random_score(np.array(p), sigma))

    return scores

def verify_model(p, sigma, N, title=None, save=False):
    model_value = expected_value.expected_value('output/triangles/points_per_mm=5', p, sigma)

    scores = get_N_random_scores(p, sigma, N)
    scores = np.array(scores)
    running_average = np.cumsum(scores) / np.arange(1, N + 1)

    if title is None:
        title = f'Monte Carlo Convergence for $\\sigma = {sigma}$ Aiming at $p = ({p[0]}, {p[1]})$'
    save_name = re.sub(' ', '', f'Monte Carlo Convergence for sigma = {sigma} Aiming at p = ({p[0]}, {p[1]})')

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, N + 1), running_average, label='Monte Carlo Running Average')
    plt.axhline(y=model_value, color='r', linestyle='--', label='Theoretical Expected Value')
    plt.xlabel('Number of Throws')
    plt.ylabel('Running Average Score')
    plt.title(title)
    plt.legend()
    plt.grid(True)

    if save:
        plt.savefig(save_name + '.png')
        plt.close()
