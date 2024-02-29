import matplotlib.pyplot as plt
import numpy as np

time = 1000
x_size = 23
y_size = 11
z_size = 101
low = -0.1
high = 0.1
y_slice = 6
x_size = 23
y_size = 11
z_size = 101

plt.rcParams["text.usetex"] = True
a = 0.19558
z_len = 1
time_step = 0.000000000001863238308179616 * 27

ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
# ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice, :]

ey_chunk_src = []
ey_chunk_out = []

for time in range(len(ey[:, 0])):
    ey_chunk_src.append(
        ey[time].reshape((x_size, y_size, z_size), order="F")[12, 6, -5]
    )
    ey_chunk_out.append(
        ey[time].reshape((x_size, y_size, z_size), order="F")[12, 6, 1]
    )


t_steps = len(ey_chunk_src)
t = np.linspace(0, t_steps - 1, t_steps)
# plt.plot(t, ey_chunk)

# plt.show()

# Compute the FFT for both signals
fft_result1 = np.fft.fft(ey_chunk_src)
fft_result2 = np.fft.fft(ey_chunk_out)
frequencies = np.fft.fftfreq(
    t_steps, d=time_step
)  # Convert time step to sampling rate

# Normalize the magnitude of signal2 by dividing it by the maximum magnitude value of signal1
fft_magnitude1 = np.abs(fft_result1)
fft_magnitude2 = np.abs(fft_result2)
max_magnitude1 = np.max(fft_magnitude1)

# Filter frequencies within the desired range (0 to 20 GHz)
desired_freq_range = (0, 15e9)  # 20 GHz
valid_indices = np.where(
    (frequencies >= desired_freq_range[0]) & (frequencies <= desired_freq_range[1])
)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(
    frequencies[valid_indices],
    fft_magnitude2[valid_indices],
    label="Outport Spectrum",
)
plt.plot(
    frequencies[valid_indices],
    fft_magnitude1[valid_indices],
    label="Source Spectrum",
)  # Add the vertical line
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Amplitude")
plt.title("Frequency Domain Representation (FFT) - Frequencies up to 20 GHz")
plt.legend()

plt.yscale("log")

plt.axvline(x=6.557e9, color="red", linestyle="--", label="Analytic Cutoff Frequency")

plt.show()
