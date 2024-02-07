from numpy import ndarray, meshgrid, linspace, mgrid, loadtxt
import plotly.graph_objects as go


def plot_volume(data: ndarray):
    # define constants
    x_size = 67
    y_size = 67
    z_size = 67
    
    # Create a figure
    fig = go.Figure()

    for t in range(data.shape[0]):
        # Reshape the 1D array into a 3D array
        # Reshape the 1D array into a 3D array
        voxels = data[t].reshape((x_size, y_size, z_size))

        # Add a volume trace for this time step
        fig.add_trace(
            go.Volume(
                x=[
                    i
                    for i in range(voxels.shape[0])
                    for _ in range(voxels.shape[1])
                    for _ in range(voxels.shape[2])
                ],
                y=[
                    j
                    for _ in range(voxels.shape[0])
                    for j in range(voxels.shape[1])
                    for _ in range(voxels.shape[2])
                ],
                z=[
                    k
                    for _ in range(voxels.shape[0])
                    for _ in range(voxels.shape[1])
                    for k in range(voxels.shape[2])
                ],
                value=voxels.flatten(),
                isomin=-1,
                isomax=1,
                opacity=0.1,  # needs to be small to see through all surfaces
                surface_count=17,  # needs to be a large number for good volume rendering
                caps=dict(x_show=False, y_show=False, z_show=False),
                surface=dict(show=True, count=2000, fill=1),
                slices=dict(x=dict(show=False), y=dict(show=False), z=dict(show=False)),
                visible=False,  # Only show this trace for its corresponding frame
            )
        )

    # Create frames for each time step
    frames = [
        go.Frame(data=[go.Volume(visible=True)], name=str(t))
        for t in range(voxels.shape[0])
    ]

    # Add frames to the figure
    fig.frames = frames

    # Create a slider to control the animation
    sliders = [
        dict(
            steps=[
                dict(
                    method="animate",
                    args=[
                        [f.name],
                        dict(frame=dict(duration=100, redraw=True), mode="immediate"),
                    ],
                    label=f.name,
                )
                for f in fig.frames
            ],
            active=0,
        )
    ]

    # Add the slider to the figure
    fig.layout.sliders = sliders

    # Show the figure
    fig.show()

def make_video(voxels_2d):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from matplotlib import cm
    from matplotlib.colors import Normalize

    # Assuming you have voxel data in a 2D numpy array
    # where each row represents a different time step

    # define constants
    size_x = 67
    size_y = 67
    size_z = 67
    vmin = -1
    vmax = 1
    cmap = cm.coolwarm

    # Create a directory to store the images
    os.makedirs('images', exist_ok=True)

    # Loop over each time step
    for t in range(voxels_2d.shape[0]):
        # Reshape the 1D array into a 3D array
        voxels = voxels_2d[t].reshape((size_x, size_y, size_z), order="F")

        # Create a 3D figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

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
        plt.savefig(f'images/frame_{t:03d}.png')

        # Close the figure to free up memory
        plt.close(fig)


def main():
    data = loadtxt("./hx.csv", delimiter=",", encoding="UTF-8")

    #plot_volume(data)
    make_video(data)


if __name__ == "__main__":
    main()
