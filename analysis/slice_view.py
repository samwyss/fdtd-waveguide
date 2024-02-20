import matplotlib.pyplot as plt
import numpy as np

time = 150
x_size = 14
y_size = 7
z_size = 67
low = -0.1
high = 0.1
slice = 4

plt.rcParams["text.usetex"] = True
a = 0.19558
z_len = 1

ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, slice, :]

fig, ax = plt.subplots()

fig.dpi = 300
fig.set_size_inches(4, 1.5)

ax.set_xlabel(r"$\hat{z}$-Position [m]")
ax.set_ylabel(r"$\hat{x}$-Position [m]")
plt.minorticks_on()

ax.imshow(ey_slice, cmap="coolwarm", vmax=high, vmin=low, extent=[0, z_len, 0, a])

plt.show()
