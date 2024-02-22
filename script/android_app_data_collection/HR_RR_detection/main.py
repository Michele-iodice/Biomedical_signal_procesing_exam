import numpy as np
from scipy.signal import butter, filtfilt, find_peaks
import csv
from scipy.interpolate import interp1d

# Assuming 'signal' is your accelerometer signal
sampling_rate = 15  # Hz
window_size = 15  # seconds


def interpolate(time, data):
    # Interpolazione lineare per ottenere dati campionati a intervalli regolari
    interpolation_function = interp1d(time, data, kind='linear', fill_value='extrapolate')
    # Definisci nuovi tempi campionati a intervalli regolari
    new_time = np.linspace(np.min(time), np.max(time), num=len(time), endpoint=True)
    # Ottieni i nuovi dati interpolati
    new_data = interpolation_function(new_time)

    return new_time, new_data


def fourier(data):
    # Calcola la trasformata di Fourier
    fft_result = np.fft.fft(data)
    fft_freq = np.fft.fftfreq(len(data), 1.0 / sampling_rate)
    return fft_freq, fft_result


def bandpass_filter(signal, lowCut, highCut, fs):
    nyquist = 0.5 * fs
    low = lowCut / nyquist
    high = highCut / nyquist
    b, a = butter(2, [low, high], btype='band')
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal


def lowpass_filter(signal, cutoff, fs):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(1, normal_cutoff, btype='low')
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal


def calculate_heart_rate(time, signal, sampling_rate):
    # Apply bandpass filter (0.5 to 3 Hz)
    interpolate_time, interpolate_signal=interpolate(time, signal)
    fft_time, fft_signal = fourier(interpolate_signal)
    filtered_signal = bandpass_filter(fft_signal, 0.5, 3, sampling_rate)
    # Find peaks in the filtered signal
    peaks, _ = find_peaks(filtered_signal)

    # Calculate heart rate in BPM
    heart_rate = round(60 / np.mean(np.diff(peaks) / sampling_rate))
    return heart_rate


def calculate_respiration_rate(signal, sampling_rate):
    # Apply low-pass filter (0.1 Hz)
    filtered_signal = lowpass_filter(signal, 0.1, sampling_rate)

    # Find peaks in the filtered signal
    peaks, _ = find_peaks(filtered_signal, height=0.5)

    # Calculate respiration rate in breaths per minute
    respiration_rate =round(60 / np.mean(np.diff(peaks) / sampling_rate))
    return respiration_rate


def load_data(filename):
    x = []
    y = []
    z = []
    time = []
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        for line in csv_reader:
            if 'X' not in line:
                elem = float(line[0])
                time.append(elem)
                elem = float(line[1])
                x.append(elem)
                elem = float(line[2])
                y.append(elem)
                elem = float(line[3])
                z.append(elem)

    return time, x, y, z


# Simulate a 15-second window of accelerometer data at 15 Hz
time = np.arange(0, window_size, 1/sampling_rate)
filename="C:/Users/39392/Desktop/University/MAGISTRALE/biomedical_signal_processing/projects/dataset/data4.csv"
time, signal_x, signal_y, signal_z = load_data(filename)

# Calculate heart rate and respiration rate for the 15-second window
heart_rate = calculate_heart_rate(time, signal_x, sampling_rate)
respiration_rate = calculate_respiration_rate(signal_z, sampling_rate)

print(f"Heart Rate: {heart_rate} BPM")
print(f"Respiration Rate: {respiration_rate} breaths/min")
