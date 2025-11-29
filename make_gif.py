import imageio.v2 as imageio

sds = ['1.0', '5.0', '10.0', '25.0', '50.0', '75.0', '100.0', '125.0', '150.0', '175.0', '200.0']
filenames = [f'StandardDev={sd}mm.png' for sd in sds]

# Target: total length = 5s
fps = len(filenames) / 5  # frames per second

with imageio.get_writer("output/normalised_heatmaps.gif", mode="I", fps=fps) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)