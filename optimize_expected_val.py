import numpy as np
import expected_value
import time
from concurrent.futures import ProcessPoolExecutor
import os
import tempfile

branch_factor = {
    0: 1,
    1: 0
}

def worker(chunk, file_to_write, triangles_dir, sd, batch_id):
    total = len(chunk)
    next_milestone = 0

    with open(file_to_write, "w") as f:
        for i, (x, y) in enumerate(chunk):
            val = expected_value.expected_value(triangles_dir=triangles_dir, p=(float(x), float(y)), std_dev=sd)
            print(f'{float(x)} {float(y)} {val}', file=f)

            progress = (i + 1) / total
            if progress >= next_milestone:
                print(f"[Worker {batch_id}] {round(next_milestone * 100)}% complete")
                f.flush()
                next_milestone += 0.5



def grid_search_subspace(sd, center, side_length, res, triangles_dir, it_layer=0):
    """
    Assumes a traingle directory is already generated
    """
    print(f'called on layer {it_layer}')
    x_range = np.linspace(center[0] - (side_length / 2), center[0] + (side_length / 2), res)
    y_range = np.linspace(center[1] - (side_length / 2), center[1] + (side_length / 2), res)

    scores = []
    X, Y = np.meshgrid(x_range, y_range)
    points = list(zip(X.ravel(), Y.ravel()))

    # separate into separate CPU's
    num_workers = os.cpu_count()
    points_chunks = [points[i::num_workers] for i in range(num_workers)]

    # Prepare filenames for temporary files
    temp_files = []
    for _ in range(num_workers):
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_files.append(temp_file.name)  # save only the name
        temp_file.close()  # close it, worker will open it

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i, chunk in enumerate(points_chunks):
            futures.append(executor.submit(worker, chunk, temp_files[i], triangles_dir, sd, i + 1))
        # Wait for all futures to complete
        for future in futures:
            future.result()

    # open the results
    for filename in temp_files:
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                x, y, exp_val = line.split(' ')
                scores.append((float(x), float(y), float(exp_val)))

    # sort the scores
    scores.sort(key=lambda tup: -tup[2]) # the expected value

    if branch_factor[it_layer] == 0:
        # just return the best tuple found
        return scores[0]
    else:
        # search around the best "couple of points" and return the best one
        local_maximas = []
        for j in range(branch_factor[it_layer]):
            local_maximas.append(grid_search_subspace(sd=sd,
                                                    center=(scores[j][0],scores[j][1]),
                                                    side_length=side_length*0.1,
                                                    res=int(res*0.25),
                                                    triangles_dir=triangles_dir,
                                                    it_layer=it_layer+1))

        local_maximas.sort(key=lambda tup: -tup[2])
        return local_maximas[0]


if __name__ == '__main__':
    # sds = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    sds = [5, 15, 25, 35, 45, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

    for sd in sds:
        print(f'-------------- Starting sd = {sd} --------------')
        start_time = time.time()
        max_point = grid_search_subspace(sd, (0,0), 230 * 2, 116, 'output/triangles/points_per_mm=2.5')
        time_taken_min = round((time.time() - start_time) / 60, 3)

        print(f'Took {time_taken_min} minutes')

        with open('max_exp_val_over_sds.csv', 'a') as f:
            print(f'{sd} {float(max_point[0])} {float(max_point[1])} {float(max_point[2])}', file=f)

        print()
