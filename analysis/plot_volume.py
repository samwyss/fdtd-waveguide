import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm
from matplotlib.colors import Normalize
from multiprocessing import Pool


def make_video(data_arrays):

    # Assuming you have voxel data in a 2D numpy array
    # where each row represents a different time step

    # define constants
    size_x = 34
    size_y = 34
    size_z = 34
    vmin = -1
    vmax = 1
    cmap = cm.coolwarm

    # Create a directory to store the images
    os.makedirs("images", exist_ok=True)

    # loop over each time step
    for t in range(data_arrays[0].shape[0]):
        # Create a 3D figure
        fig = plt.figure(figsize=(10, 10))

        # Loop over each spatial index
        for i, voxels_2d in enumerate(data_arrays):
            # Reshape the 1D array into a 3D array
            voxels = voxels_2d[t].reshape((size_x, size_y, size_z), order="F")

            # Create a 3D subplot
            ax = fig.add_subplot(2, 3, i + 1, projection="3d")

            # Get the coordinates of the voxels
            x, y, z = np.indices((size_x + 1, size_y + 1, size_z + 1))

            # Create a mask that is True for all voxels
            mask = np.ones((size_x, size_y, size_z), dtype=bool)

            # Create a color array using a colormap
            colors = cmap((voxels - vmin) / (vmax - vmin))

            # Draw the voxels
            ax.voxels(x, y, z, mask, facecolors=colors)

            # Create a ScalarMappable object for the color bar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
            sm.set_array([])

            # Add the color bar to the figure
            fig.colorbar(sm, ax=ax)

        # Save the figure as an image
        plt.savefig(f"images/frame_{t:03d}.png")

        # Close the figure to free up memory
        plt.close(fig)


def make_frames(data_arrays):
    # Function to render a single frame

    # Assuming you have six 2D numpy arrays
    # where each row represents a different time step

    # Create a directory to store the images
    os.makedirs("images", exist_ok=True)

    # Create a pool of worker processes
    with Pool() as pool:
        # Use the pool to render the frames in parallel
        pool.map(
            render_frame, [(t, data_arrays) for t in range(data_arrays[0].shape[0])]
        )


def main():
    hx = np.loadtxt("./hx.csv", delimiter=",", encoding="UTF-8")
    hy = np.loadtxt("./hy.csv", delimiter=",", encoding="UTF-8")
    hz = np.loadtxt("./hz.csv", delimiter=",", encoding="UTF-8")
    ex = np.loadtxt("./ex.csv", delimiter=",", encoding="UTF-8")
    ey = np.loadtxt("./ey.csv", delimiter=",", encoding="UTF-8")
    ez = np.loadtxt("./ez.csv", delimiter=",", encoding="UTF-8")

    data = [hx, hy, hz, ex, ey, ez]

    print("loaded data rendering now")

    # plot_volume(data)
    # make_video(data)
    make_frames(data)


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from multiprocessing import Pool


# Function to render a single frame
def render_frame(args):
    # Set the fixed range for the colormap
    size_x = 17
    size_y = 17
    size_z = 17
    vmine = -0.01
    vmaxe = 0.01
    vminh = -0.00001
    vmaxh = 0.00001
    cmap = cm.coolwarm
    t, data_arrays = args
    # Create a 3D figure
    fig = plt.figure(figsize=(10, 10))

    for i, voxels_2d in enumerate(data_arrays):

        # Create a 3D subplot
        ax = fig.add_subplot(2, 3, i + 1, projection="3d")
        ax.view_init(45, 45)

        ax.set_xlabel("X Axis")
        ax.set_ylabel("Y Axis")
        ax.set_zlabel("Z Axis")

        if i == 0:
            ax.set_title("Hx")
        elif i == 1:
            ax.set_title("Hy")
        elif i == 2:
            ax.set_title("Hz")
        elif i == 3:
            ax.set_title("Ex")
        elif i == 4:
            ax.set_title("Ey")
        elif i == 5:
            ax.set_title("Ez")

        # Reshape the 1D array into a 3D array
        voxels = voxels_2d[t].reshape((size_x, size_y, size_z), order="F")

        # Get the coordinates of the voxels
        x, y, z = np.indices((size_x + 1, size_y + 1, size_z + 1))

        # Create a mask that is True for all voxels
        mask = np.ones((size_x, size_y, size_z), dtype=bool)

        # Create a color array using a colormap

        if i >= 2 == 1:
            colors = cmap((voxels - vmine) / (vmaxe - vmine))
        else:
            colors = cmap((voxels - vminh) / (vmaxh - vminh))

        # Draw the voxels
        ax.voxels(x, y, z, mask, facecolors=colors)

        # Create a ScalarMappable object for the color bar
        if i >= 3:
            sm = plt.cm.ScalarMappable(
                cmap=cmap, norm=Normalize(vmin=vmine, vmax=vmaxe)
            )
        else:
            sm = plt.cm.ScalarMappable(
                cmap=cmap, norm=Normalize(vmin=vminh, vmax=vmaxh)
            )
        sm.set_array([])

        # Add the color bar to the figure
        # fig.colorbar(sm, ax=ax)

    # Save the figure as an image
    plt.savefig(f"images/frame_{t:03d}.png")

    # Close the figure to free up memory
    plt.close(fig)


if __name__ == "__main__":
    main()
