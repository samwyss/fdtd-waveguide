from matplotlib import streamplot
import matplotlib.pyplot as plt
import numpy as np

time = 1000 
x_size = 23
y_size = 11
z_size = 101
elow = -0.1
ehigh = 0.1
hlow = -0.0001
hhigh = 0.0001
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
hy_slice2 = hy[time].reshape((x_size, y_size, z_size), order="F")[x_slice_idx, :, :]

x = np.linspace(0, a, x_size)
y = np.linspace(0,b, y_size)
z = np.linspace(0, z_len, z_size)

fig, ax = plt.subplots(2,1)
fig.dpi = 300
fig.set_size_inches(4, 1.5)

""" ax.set_xlabel(r"$\hat{z}$-Position [m]")
ax.set_ylabel(r"$\hat{x}$-Position [m]") """
plt.minorticks_on()

X1, Z1 = np.meshgrid(z, x)
Y2, Z2 = np.meshgrid(z, y)


ax[0].quiver(
    X1, Z1, hz_slice, hx_slice, scale=0.04, pivot="tip", angles="uv"
)  # Adjust the scale as needed
ax[0].imshow(ey_slice, cmap="coolwarm", vmax=ehigh, vmin=elow, extent=[0, z_len, 0, a])
# ax.streamplot(X, Z, hz_slice, hx_slice, density=1, color='k')

ax[1].quiver(
Y2, Z2, ez_slice, ey_slice2, scale=10, pivot="tip", angles="uv"
)  # Adjust the scale as needed
ax[1].imshow(hy_slice2, cmap="coolwarm", vmax=hhigh, vmin=hlow, extent=[0, z_len, 0, b])
# ax.streamplot(X, Z, hz_slice, hx_slice, density=1, color='k')
plt.show()
