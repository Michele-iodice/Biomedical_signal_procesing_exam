import numpy as np
from scipy.signal import find_peaks, butter, filtfilt


def extract_acceleration_magnitude(data):
    # Extract the magnitude of acceleration from X, Y, Z axes
    return np.sqrt(np.sum(np.square(data[:, 1:]), axis=1))


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


def calculate_heart_rate(timestamp_acceleration_data, sampling_rate):
    # Extract acceleration signal
    acceleration_data = np.array(timestamp_acceleration_data)
    acceleration_magnitude = extract_acceleration_magnitude(acceleration_data)

    # Apply bandpass filter (0.5 to 3 Hz)
    filtered_signal = bandpass_filter(acceleration_magnitude, 0.5, 3, sampling_rate)

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


# Similarly, you can modify the function for respiration rate
# ...

# Example usage:
# Assume 'timestamp_acceleration_data' is your timestamped accelerometer data
sampling_rate = 15
average_heart_rate = calculate_heart_rate(timestamp_acceleration_data, sampling_rate)

print("Average Heart Rate (BPM):", average_heart_rate)
