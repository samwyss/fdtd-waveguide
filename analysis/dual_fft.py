import matplotlib.pyplot as plt
import numpy as np

time = 50
x_size = 16
y_size = 7
z_size = 34
low = -0.1
high = 0.1
y_slice = 6

plt.rcParams["text.usetex"] = True
a = 0.19558
z_len = 1
time_step = 0.000000000002792360102758852

eybee = np.loadtxt("eybee.csv", delimiter=",", encoding="UTF-8")
eyber = np.loadtxt("eyber.csv", delimiter=",", encoding="UTF-8")


ey_chunk_1 = []
ey_chunk_2 = []
""" ey_chunk_ref = [] """

for time in range(len(eybee[:, 0])):
    ey_chunk_1.append(
        eybee[time].reshape((29, 13, 62), order="F")[15, 7, 31]
    )

for time in range(len(eyber[:, 0])):
    ey_chunk_2.append(
        eyber[time].reshape((47, 21, 102), order="F")[24, 11, 51]
    )


t_steps1 = len(ey_chunk_1)
t_steps2 = len(ey_chunk_2)
t1 = np.linspace(0, t_steps1 - 1, t_steps1)
t2 = np.linspace(0, t_steps2-1, t_steps2)
# plt.plot(t, ey_chunk)

# plt.show()

# Compute the FFT for both signals
fft_result1 = np.fft.fft(ey_chunk_1)
fft_result2 = np.fft.fft(ey_chunk_2)

frequencies1 = np.fft.fftfreq(t_steps1, d=6.1e-12)
frequencies2 = np.fft.fftfreq(t_steps2, d=5.6241e-12)

# Normalize the magnitude of signal2 by dividing it by the maximum magnitude value of signal1
fft_magnitude1 = np.abs(fft_result1)
fft_magnitude2 = np.abs(fft_result2)

# Filter frequencies within the desired range (0 to 20 GHz)
desired_freq_range = (1, 12e9)  # 20 GHz
valid_indices1 = np.where(
    (frequencies1 >= desired_freq_range[0]) & (frequencies1 <= desired_freq_range[1])
)
valid_indices2 = np.where(
    (frequencies2 >= desired_freq_range[0]) & (frequencies2 <= desired_freq_range[1])
)

plt.rcParams["figure.dpi"] = 300

# Plot the results
plt.figure()
plt.plot(
    frequencies1[valid_indices1],
    fft_magnitude2[valid_indices1],
    label="Beeswax",
)
plt.plot(
    frequencies2[valid_indices2],
    fft_magnitude1[valid_indices2],
    label="Beryllia",
)

plt.xlabel("Frequency (Hz)")
plt.ylabel("Intenisty (Arb. Units)")
plt.legend()

plt.yscale("log")
plt.xlim((1e9, 12e9))

# plt.axvline(x=6.557e9, color="red", linestyle="--", label="Analytic Cutoff Frequency")

plt.minorticks_on()
plt.tick_params(
    which="both",
    axis="both",
    top=True,
    right=True,
    labeltop=False,
    labelright=False,
)

plt.legend()

plt.show()
