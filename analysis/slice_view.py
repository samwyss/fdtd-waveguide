import matplotlib.pyplot as plt
import numpy as np

time = 50
x_size = 33
y_size = 17
z_size = 334
low = -0.1
high = 0.1
slice = 9

plt.rcParams["text.usetex"] = True

ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, 9, :]

fig, ax = plt.subplots()

fig.dpi = 300
fig.set_size_inches(4, 1.5)

ax.set_xlabel(r"$\hat{z}$-Position [m]")
ax.set_ylabel(r"$\hat{x}$-Position [m]")
plt.minorticks_on()

ax.imshow(ey_slice, cmap="coolwarm", vmax=high, vmin=low, extent=[0, 2, 0, 0.19558])

plt.show()
