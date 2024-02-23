import matplotlib.pyplot as plt
import numpy as np

time = 150
x_size = 14
y_size = 7
z_size = 67
low = -0.1
high = 0.1
y_slice = 4
x_size = 23
y_size = 11
z_size = 101

plt.rcParams["text.usetex"] = True
a = 0.19558
z_len = 1
time_step = 0.000000000001863238308179616

ey = np.loadtxt("./out/ey.csv", delimiter=",", encoding="UTF-8")
# ey_slice = ey[time].reshape((x_size, y_size, z_size), order="F")[:, y_slice, :]

ey_chunk = []

for time in range(len(ey[:, 0])):
    ey_chunk.append(ey[time].reshape((x_size, y_size, z_size), order="F")[12, 6, 1])

t_steps = len(ey_chunk)
t = np.linspace(0, t_steps - 1, t_steps)
# plt.plot(t, ey_chunk)

# plt.show()

# Compute the FFT
fft_result = np.fft.fft(ey_chunk)
frequencies = np.fft.fftfreq(
    t_steps, d=time_step
)  # Convert time step to sampling rate

# Plot the results
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(t, ey_chunk)
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Original Time Series Data")

plt.subplot(2, 1, 2)
plt.plot(frequencies, np.abs(fft_result))
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.title("Frequency Domain Representation (FFT)")

plt.tight_layout()
plt.show()

# Compute the FFT
fft_result = np.fft.fft(ey_chunk)
frequencies = np.fft.fftfreq(
    t_steps, d=time_step
)  # Convert time step to sampling rate

# Filter frequencies within the desired range (0 to 20 GHz)
desired_freq_range = (0, 15e9)  # 20 GHz
valid_indices = np.where(
    (frequencies >= desired_freq_range[0]) & (frequencies <= desired_freq_range[1])
)

# Normalize the time series data
signal_normalized = (ey_chunk - np.mean(ey_chunk)) / np.std(ey_chunk)

# Normalize the FFT result (magnitude)
fft_magnitude = np.abs(fft_result)
fft_magnitude_normalized = fft_magnitude / np.max(fft_magnitude)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(frequencies[valid_indices], fft_magnitude_normalized[valid_indices])
plt.xlabel("Frequency (Hz)")
plt.ylabel("Normalized Amplitude")
plt.title("Frequency Domain Representation (FFT) - Frequencies up to 20 GHz")
plt.xscale("log")

plt.tight_layout()
plt.grid(True)
plt.show()