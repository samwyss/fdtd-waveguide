using CairoMakie
using DelimitedFiles

# load data
data = readdlm("./out/ey.csv", ',', Float64, '\n')

# constants
time = 50
x_size = 33
y_size = 17
z_size = 334
low = -0.1
high = 0.1

# slice data on timestep
data_slice = data[time, :]

# reorder into volume
data_vol = reshape(data_slice, (x_size, y_size, z_size))

f = Figure()
ax = LScene(f[1, 1], show_axis=false)
Axis(f[1, 1])
ax.aspect = x_size / z_size

co = contourf!(data_vol[:, 9, :], levels = range(low, high, length = 6), colormap=:coolwarm, colorrange=(low,high))#interpolate=true

Colorbar(f[1, 2], co, ticks = low:0.2:high)

#colsize!(f.layout, 1, Aspect(1, x_size / z_size))

#resize_to_layout!(f)

f