from matplotlib import streamplot
import matplotlib.pyplot as plt
import numpy as np

time = 628 
x_size = 23
y_size = 11
z_size = 101
elow = -0.1
ehigh = 0.1
hlow = -0.000351
hhigh = 0.000351
y_slice_idx = 6
x_slice_idx = 12

plt.rcParams["text.usetex"] = True
a = 0.02286
b=0.01016
z_len = 0.1

ex = np.loadtxt("./out/ex.csv", delimiter=",", encoding="UTF-8")
ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
ez = np.loadtxt("./out/ez.csv", delimiter=",", encoding="UTF-8")
hx = np.loadtxt("./out/hx.csv", delimiter=",", encoding="UTF-8")
hy = np.loadtxt("./out/hy.csv", delimiter=",", encoding="UTF-8")
hz = np.loadtxt("./out/hz.csv", delimiter=",", encoding="UTF-8")

ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice_idx, :]
hx_slice = hx[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice_idx, :]
hz_slice = hz[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice_idx, :]

ey_slice2 = ey[time].reshape((x_size, y_size, z_size), order="F")[x_slice_idx, :, :]
ez_slice = ez[time].reshape((x_size, y_size, z_size), order="F")[x_slice_idx, :, :]
hx_slice2 = hx[time].reshape((x_size, y_size, z_size), order="F")[x_slice_idx, :, :]

x = np.linspace(0, a, x_size)
y = np.linspace(0,b, y_size)
z = np.linspace(0, z_len, z_size)


plt.rcParams["figure.dpi"] = 600


fig, ax = plt.subplots(2,1)

ax[0].set_ylabel(r"$a$ [m]")
plt.minorticks_on()

X1, Z1 = np.meshgrid(z, x)
Y2, Z2 = np.meshgrid(z, y)


ax[0].quiver(
    X1, Z1, hz_slice, hx_slice, scale=0.04, pivot="mid", angles="uv"
)  # Adjust the scale as needed
im0 = ax[0].imshow(ey_slice, cmap="coolwarm", vmax=ehigh, vmin=elow, extent=[0, z_len, 0, a])
# ax.streamplot(X, Z, hz_slice, hx_slice, density=1, color='k')
""" cbar0 = plt.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
cbar0.set_label("Intensity [V/m]") """
# ax[0].ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
ax[0].minorticks_on()
ax[0].tick_params(
    which="both",
    axis="both",
    top=True,
    right=True,
    labeltop=False,
    labelright=False,
)

ax[1].set_xlabel(r"$l$ [m]")
ax[1].set_ylabel(r"$b$ [m]")

ax[1].quiver(
Y2, Z2, ez_slice, ey_slice2, scale=15, pivot="mid", angles="uv"
)  # Adjust the scale as needed
im1 = ax[1].imshow(hx_slice2, cmap="coolwarm", vmax=hhigh, vmin=hlow, extent=[0, z_len, 0, b])
# ax.streamplot(X, Z, hz_slice, hx_slice, density=1, color='k')
""" cbar0 = plt.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
cbar0.set_label("Intensity [A/m]") """
# ax[1].ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
ax[1].minorticks_on()
ax[1].tick_params(
    which="both",
    axis="both",
    top=True,
    right=True,
    labeltop=False,
    labelright=False,
)

plt.show()
