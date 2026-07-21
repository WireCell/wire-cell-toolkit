#!/usr/bin/env -S uv run --script
# -*- python -*-
# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "matplotlib","scipy","pyqt5"]
# ///
'''
This arose from a fight between me and Gemini due to an inability to make
clear I was asking about how properly roll a convolution result and assign a
correct time to the first sample of the rolled result.  

If the first sample of the measure array is t=0 then we want the first sample of
the final convolution to have negative time.  This is because any "delay" in the
deconvolved response corresponds to a "backing up in time" of the signal.  

In any case, I'm right, Gemini is wrong, so there!

If we roll by N we label the rolled sample zero by t = t0 - dt*N where t0 is the
time of the zero sample of the measure.

Since we pad by Nr-1 for the response and Nf-1 for the filter to give a total
size of Nm+Nr+Nf-2, I believe we should roll by Nr-1.
'''



import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import general_gaussian

# --- Parameters ---
dt = 1.0  # Time step (e.g., 1 unit). Represents the duration of one column index.

# Lengths of original 1D signals (representing columns of 2D arrays)
len_m = 200 # Length of measure m
len_f = 20 # Length of filter f
#len_r = 25 # Length of detector response r
len_r = 100 # Length of detector response r

# Define padding length as specified: sum of lengths - 2
# This length is generally suitable for linear deconvolution via FFT.
pad_len = len_m + len_f + len_r - 2

# --- Signal Generation (1D for simplicity, represents one column from a 2D array) ---
# We assume 'm[0]' corresponds to physical time t=0 for simplicity in labeling.

# m: Measure with a clear peak
m_peak_idx = 10 # Peak at index 10 (i.e., physical time 10*dt)
m = np.zeros(len_m)
sigma_m = 2.0
m[m_peak_idx - int(2*sigma_m) : m_peak_idx + int(2*sigma_m) + 1] = \
    np.exp(-0.5 * ((np.arange(-int(2*sigma_m), int(2*sigma_m) + 1)) / sigma_m)**2)

# f: Filter with a delayed peak (introduces a positive shift in the deconvolution result)
f_peak_idx = 15 # Peak at index 15 (i.e., effective delay of 15*dt)
f = np.zeros(len_f)
sigma_f = 2.0
f[f_peak_idx - int(2*sigma_f) : f_peak_idx + int(2*sigma_f) + 1] = \
    np.exp(-0.5 * ((np.arange(-int(2*sigma_f), int(2*sigma_f) + 1)) / sigma_f)**2)


# r: Detector response with a delayed peak (introduces a negative shift in deconvolution)
# r_peak_idx = 40 # Peak at index 5 (i.e., effective delay of 5*dt)
r_peak_idx = len_r - 10
r = np.zeros(len_r)
sigma_r = 2.0
r[r_peak_idx - int(2*sigma_r) : r_peak_idx + int(2*sigma_r) + 1] = \
    np.exp(-0.5 * ((np.arange(-int(2*sigma_r), int(2*sigma_r) + 1)) / sigma_r)**2)


# --- Deconvolution Function ---
def filtered_deconvolve(m_in, f_in, r_in, pad_len):
    # 2. Zero-pad m, f, r to the calculated size
    m_padded = np.pad(m_in, (0, pad_len - len(m_in)), 'constant')
    f_padded = np.pad(f_in, (0, pad_len - len(f_in)), 'constant')
    r_padded = np.pad(r_in, (0, pad_len - len(r_in)), 'constant')

    # 3. Apply 1D discrete Fourier transform on each
    M_freq = np.fft.fft(m_padded)
    F_freq = np.fft.fft(f_padded)
    R_freq = np.fft.fft(r_padded)

    # 4. Perform element-wise deconvolution S = M*F/R
    S_freq = np.zeros_like(M_freq, dtype=complex)
    
    # Handle R_freq elements that are zero or very small.
    # For a robust deconvolution, consider a Wiener filter or adding a small epsilon.
    # Here, we set S_freq to zero where R_freq is effectively zero to avoid division by zero.
    threshold = 1e-6 * np.abs(R_freq).max() # A small threshold relative to max R_freq
    
    non_zero_R_indices = np.abs(R_freq) > threshold
    S_freq[non_zero_R_indices] = (M_freq[non_zero_R_indices] * F_freq[non_zero_R_indices]) / R_freq[non_zero_R_indices]
    # S_freq[~non_zero_R_indices] will remain zero as initialized

    # 5. Apply inverse DFT and take the real part
    s_raw = np.fft.ifft(S_freq).real

    return s_raw

# --- Perform Deconvolution ---
s_raw = filtered_deconvolve(m, f, r, pad_len)

# --- Analysis of Shifts ---
# Find peak indices (effective delays) for m, f, r within their original lengths.
# np.argmax finds the index of the maximum value.
P_m = np.argmax(m)
P_f = np.argmax(f)
P_r = np.argmax(r)

# Calculate the expected shift of m's peak in s (in samples)
# Net shift = (Delay from f) - (Delay from r)
expected_net_shift_samples = P_f - P_r
expected_peak_index_in_s = P_m + expected_net_shift_samples

# Find the observed peak index in s_raw (within the padded output)
# We look within a reasonable range to avoid noise/ringing peaks far from expected.
# Using the full range of s_raw:
P_s_raw = np.argmax(s_raw) 
observed_shift_samples = P_s_raw - P_m

print(f"--- Shift Analysis ---")
print(f"Original m peak at index (P_m): {P_m}")
print(f"Filter f peak at index (P_f): {P_f}")
print(f"Detector r peak at index (P_r): {P_r}")
print(f"Expected net shift (P_f - P_r): {expected_net_shift_samples} samples")
print(f"Expected s peak index (P_m + net shift): {expected_peak_index_in_s}")
print(f"Observed s_raw peak at index: {P_s_raw}")
print(f"Observed shift from m's peak: {observed_shift_samples}")
print(f"Difference (Expected - Observed): {expected_peak_index_in_s - P_s_raw} (should be ~0 for ideal conditions)")

# Calculate physical time of s_raw[0] assuming m[0] is at physical time T_m0=0.
# The deconvolution effectively maps the original 't=0' of m to 't=-(P_f - P_r)' in the output's time frame.
t_physical_s_raw_0 = -(P_f - P_r) * dt
print(f"Physical time represented by s_raw[0] (assuming m[0] is at t=0): {t_physical_s_raw_0:.2f} s")

# --- Plotting ---
# Define time axes for plotting
time_axis_m = np.arange(len_m) * dt
time_axis_f = np.arange(len_f) * dt
time_axis_r = np.arange(len_r) * dt
time_axis_s_raw = np.arange(pad_len) * dt # The x-axis for s_raw, where s_raw[k] is at k*dt

fig, axs = plt.subplots(3, 1, figsize=(12, 12))

# Plot m, f, r
axs[0].plot(time_axis_m, m, label='m (Measure)')
axs[0].plot(time_axis_f, f, label='f (Filter)')
axs[0].plot(time_axis_r, r, label='r (Detector Response)')
axs[0].axvline(P_m * dt, color='blue', linestyle='--', label=f'm peak ({P_m*dt:.1f}s)')
axs[0].axvline(P_f * dt, color='orange', linestyle='--', label=f'f peak ({P_f*dt:.1f}s)')
axs[0].axvline(P_r * dt, color='green', linestyle='--', label=f'r peak ({P_r*dt:.1f}s)')
axs[0].set_title('1. Input Signals (m, f, r) and their Peaks')
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Amplitude')
axs[0].legend()
axs[0].grid(True)

# Plot s_raw and compare with m
axs[1].plot(time_axis_m, m, label='m (Measure)', linestyle=':', color='gray', alpha=0.7)
axs[1].plot(time_axis_s_raw, s_raw, label='s (Deconvolved Output)', color='red')
axs[1].axvline(P_m * dt, color='gray', linestyle='--', label=f'm peak ({P_m*dt:.1f}s)')
axs[1].axvline(P_s_raw * dt, color='red', linestyle='--', label=f's peak ({P_s_raw*dt:.1f}s)')
axs[1].axvline(expected_peak_index_in_s * dt, color='purple', linestyle=':', label=f'Expected s peak ({expected_peak_index_in_s*dt:.1f}s)', alpha=0.8)
axs[1].set_title(f'2. Deconvolution Result "s" vs. "m" (Observed shift: {observed_shift_samples}, Expected: {expected_net_shift_samples})')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Amplitude')
axs[1].legend()
axs[1].grid(True)
axs[1].set_xlim(0, max(time_axis_m[-1], time_axis_s_raw[P_s_raw] + 20*dt)) # Adjust xlim to clearly see shifts

# --- Demonstrate User's Proposed Rolling ---
# User's proposed roll: cycle Nr (len_r) of the highest columns of s to the front.
# In numpy, circshift(array, -N) moves elements from index N onwards to the beginning.
# This means s_rolled[0] will contain s_raw[len_r].
#s_rolled_user = np.roll(s_raw, -len_r)

N_roll = len_r
print(f'ROLLING by {N_roll}')
s_rolled_user = np.roll(s_raw, N_roll) 

t_m_0 = 0
t_u_0 = t_m_0 - N_roll*dt

# What is the actual physical time of s_rolled_user[0]?
# s_rolled_user[0] is the value that was originally at s_raw[len_r].
# The physical time for s_raw[k] is t_physical_s_raw_0 + k*dt.
# So, for s_rolled_user[0] (which came from s_raw[len_r]), its physical time is:
t_physical_s_rolled_user_0 = t_physical_s_raw_0 + len_r * dt

# User's proposed label for s_rolled_user[0]: T - Tr.
# Assuming T=0 (m[0] is at physical time 0).
# Tr is the duration of r, which for an `len_r` length array is (len_r - 1) * dt.
time_axis_s_user = np.arange(pad_len) * dt + t_u_0

# peak index and peak time
P_u = np.argmax(s_rolled_user)
t_u_P = t_u_0 + P_u*dt


print(f"User's proposed roll: circular shift s_raw by {N_roll}")
print(f"Rolled start time: {t_u_0:.2f} s")
print(f"Rolled peak time: {t_u_P:.2f} s")




# Plot s_rolled_user
axs[2].plot(time_axis_m, m, label='m (Measure)', linestyle=':', color='gray', alpha=0.7)
axs[2].plot(time_axis_s_user, s_rolled_user, label='s_rolled (User Proposed Roll)', color='blue')
# Add a vertical line for the true physical time of s_rolled_user[0]
axs[2].axvline(t_u_P, color='purple', linestyle=':', label=f's_rolled[0] time ({t_u_P:.2f}s)')
# Add a vertical line for the user's proposed label
axs[2].set_title(f'3. "s" with User Proposed Roll (-len_r) and Labeling')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Amplitude')
axs[2].legend()
axs[2].grid(True)
axs[2].set_xlim(t_u_0, max(time_axis_m[-1], time_axis_s_raw[P_s_raw] + 20*dt))

plt.tight_layout()
plt.show()
