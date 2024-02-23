using GLMakie
using DelimitedFiles

# load data
data = readdlm("./out/ey.csv", ',', Float64, '\n')

# constants
#time = 4
x_size = 23
y_size = 11
z_size = 101

fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)

x = LinRange(1, x_size, x_size)
y = LinRange(1, y_size, y_size)
z = LinRange(1, z_size, z_size)

time = Observable(1)

sgrid = SliderGrid(
    fig[2, 1],
    (label = "yz plane - x axis", range = 1:length(x)),
    (label = "xz plane - y axis", range = 1:length(y)),
    (label = "xy plane - z axis", range = 1:length(z)),
    (label = "Time", range = 1:length(data[:,1])),
)

lo = sgrid.layout
nc = ncols(lo)

# slice data on timestep
data_slice = @lift(data[$time, :])

# reorder into volume
data_vol = @lift(reshape($data_slice, (x_size, y_size, z_size)))

low = -0.1
high = 0.1

vol = data_vol#[cos(X)*sin(Y)*sin(Z) for X ∈ x, Y ∈ y, Z ∈ z]
plt = volumeslices!(ax, x, y, z, vol, colormap=:coolwarm, colorrange=(low,high)) #interpolate=true

# connect sliders to `volumeslices` update methods
sl_yz, sl_xz, sl_xy = sgrid.sliders

on(sl_yz.value) do v; plt[:update_yz][](v) end
on(sl_xz.value) do v; plt[:update_xz][](v) end
on(sl_xy.value) do v; plt[:update_xy][](v) end
on(sgrid.sliders[4].value) do v 
    time[] = v
    # slice data on timestep
    data_slice = @lift(data[$time, :])

    # reorder into volume
    data_vol = @lift(reshape($data_slice, (x_size, y_size, z_size)))

    #vol = data_vol#[cos(X)*sin(Y)*sin(Z) for X ∈ x, Y ∈ y, Z ∈ z]
    #plt = volumeslices!(ax, x, y, z, vol, colormap=:coolwarm, colorrange=(-0.1,0.1)) # interpolate=true

    set_close_to!(sl_yz, .5length(x))
    set_close_to!(sl_xz, .5length(y))
    set_close_to!(sl_xy, .5length(z))

    # add toggles to show/hide heatmaps
    hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
    toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

    map(zip(hmaps, toggles)) do (h, t)
        connect!(h.visible, t.active)
    end
end

set_close_to!(sl_yz, .5length(x))
set_close_to!(sl_xz, .5length(y))
set_close_to!(sl_xy, .5length(z))

# add toggles to show/hide heatmaps
hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

map(zip(hmaps, toggles)) do (h, t)
    connect!(h.visible, t.active)
end

cam3d!(ax.scene, projectiontype=Makie.Orthographic)

# cam3d!(ax.scene, projectiontype=Makie.Orthographic)

fig