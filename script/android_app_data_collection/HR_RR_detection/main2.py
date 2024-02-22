import numpy as np
from scipy.signal import find_peaks, butter, filtfilt


def bandpass_filter(signal, lowcut, highcut, fs):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(1, [low, high], btype='band')
    return filtfilt(b, a, signal)


def lowpass_filter(signal, cutoff, fs):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(1, normal_cutoff, btype='low')
    return filtfilt(b, a, signal)


def calculate_heart_rate(signal, sampling_rate):
    # Apply bandpass filter (0.5 to 3 Hz)
    filtered_signal = bandpass_filter(signal, 0.5, 3, sampling_rate)

    # Calculate heart rate in BPM for 15-second windows
    window_size = 15  # seconds
    window_samples = int(window_size * sampling_rate)

    heart_rates = []
    for i in range(0, len(filtered_signal), window_samples):
        window = filtered_signal[i:i+window_samples]

        # Find peaks in the windowed signal
        peaks, _ = find_peaks(window, height=0.5)

        # Calculate heart rate in BPM
        if len(peaks) >= 2:
            heart_rate = 60 / np.mean(np.diff(peaks) / sampling_rate)
            heart_rates.append(heart_rate)

    # Return the average heart rate
    return np.mean(heart_rates)


def calculate_respiration_rate(signal, sampling_rate):
    # Apply low-pass filter (0.1 to 0.5 Hz)
    filtered_signal = lowpass_filter(signal, 0.5, sampling_rate)

    # Calculate respiration rate in breaths per minute for 15-second windows
    window_size = 15  # seconds
    window_samples = int(window_size * sampling_rate)

    respiration_rates = []
    for i in range(0, len(filtered_signal), window_samples):
        window = filtered_signal[i:i+window_samples]

        # Find peaks in the windowed signal
        peaks, _ = find_peaks(window, height=0.5)

        # Calculate respiration rate in breaths per minute
        if len(peaks) >= 2:
            respiration_rate = 60 / np.mean(np.diff(peaks) / sampling_rate)
            respiration_rates.append(respiration_rate)

    # Return the average respiration rate
    return np.mean(respiration_rates)


# Example usage:
# Assume 'acceleration_data' is your accelerometer signal with a sampling rate of 15 Hz
sampling_rate = 15
average_heart_rate = calculate_heart_rate(acceleration_data, sampling_rate)
average_respiration_rate = calculate_respiration_rate(acceleration_data, sampling_rate)

print("Average Heart Rate (BPM):", average_heart_rate)
print("Average Respiration Rate (BPM):", average_respiration_rate)
