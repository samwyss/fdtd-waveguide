from matplotlib import streamplot
import matplotlib.pyplot as plt
import numpy as np

time = 120
x_size = 23
y_size = 11
z_size = 101
low = -0.0001
high = 0.0001
y_slice = 6

plt.rcParams["text.usetex"] = True
a = 0.19558
z_len = 1

ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
hx = np.loadtxt("./out/hx.csv", delimiter=",", encoding="UTF-8")
hz = np.loadtxt("./out/hz.csv", delimiter=",", encoding="UTF-8")
ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice, :]
hx_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice, :]
hz_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice, :]

x = np.linspace(0, a, x_size)
z = np.linspace(0, z_len, z_size)

fig, ax = plt.subplots()
fig.dpi = 300
fig.set_size_inches(4, 1.5)

ax.set_xlabel(r"$\hat{z}$-Position [m]")
ax.set_ylabel(r"$\hat{x}$-Position [m]")
plt.minorticks_on()

x = np.linspace(0, a, x_size)
z = np.linspace(0, z_len, z_size)
X, Z = np.meshgrid(z, x)
plt.quiver(
    X, Z, hx_slice, hz_slice, scale=5, pivot="tip", angles="uv"
)  # Adjust the scale as needed
# ax.imshow(ey_slice, cmap="coolwarm", vmax=high, vmin=low, extent=[0, z_len, 0, a])
plt.show()

print(np.allclose(hx_slice, hz_slice))
