% Load the CSV file
file_path = 'C:/Users/39392/Documents/GitHub/Biomedical_signal_procesing_exam/dataset/data4.csv';

% Load data 
data_table = readtable(file_path);

% estraggo i dati sui tre assi x, y, e z
acc_x = data_table.X;
acc_y = data_table.Y;
acc_z = data_table.Z;

% timestamps
timestamps = data_table.TIME; % Assuming the column name is 'timestamp'

% Parametri
fs = 1 / (mean(diff(timestamps))); % Frequenza di campionamento (Hz)
duration = 16; % durata del segnale in secondi
n = length(acc_z); %lunghezza del segnale in numero di campioni
f_cutoff_resp = 2; % Frequenza di taglio del filtro passa-basso per il tasso respiratorio (Hz)
f_cutoff_heart_low = 0.1; % Frequenza di taglio inferiore del filtro passa-banda per il tasso cardiaco (Hz)
f_cutoff_heart_high = 10; % Frequenza di taglio superiore del filtro passa-banda per il tasso cardiaco (Hz)

% Filtraggio dei dati
[b_low, a_low] = butter(2, f_cutoff_resp / (fs/2), 'low'); % Coefficienti del filtro passa-basso
resp_filtered = filtfilt(b_low, a_low, acc_z); % Applicazione del filtro passa-basso

[b_band, a_band] = butter(2, [f_cutoff_heart_low, f_cutoff_heart_high] / (fs/2), 'bandpass'); % Coefficienti del filtro passa-banda
heart_filtered = filtfilt(b_band, a_band, acc_x); % Applicazione del filtro passa-banda


% Calcola la Trasformata di Fourier
frequencies = linspace(0, fs/2, n/2+1);
resp_fourier_transform = fft(resp_filtered) / n;
resp_amplitude_spectrum = abs(resp_fourier_transform(1:n/2+1));
heart_fourier_transform = fft(heart_filtered) / n;
heart_amplitude_spectrum = abs(heart_fourier_transform(1:n/2+1));


% Individuazione dei picchi respiratori
[resp_peaks, locs_peaks] = findpeaks(resp_amplitude_spectrum);



% Individuazione dei picchi cardiaci
[heart_peaks, locs_peaks_hearth] = findpeaks(heart_amplitude_spectrum);

% Calcola l'hearth rate in Hz
resp_periods = length(resp_peaks) / duration;
resp_rate = resp_periods * 60/duration;


% Calcola l'hearth rate in Hz
heart_period = length(heart_peaks) / duration;
heart_rate = heart_period * 60/duration;


disp(['Respiration rate: ', num2str(round(resp_rate)), ' bpm']);
disp(['Hear rate: ', num2str(round(heart_rate)), ' bpm']);

figure;
subplot(2, 1, 1);
plot(timestamps, acc_x, 'r');
title('Accelerometer Data along X axes');
xlabel('Time');
ylabel('X');

subplot(2, 1, 2);
plot(timestamps, acc_z, 'b');
title('Accelerometer Data along  Z axes');
xlabel('Time');
ylabel('Z');


figure;
subplot(2, 1, 1);
plot(timestamps, resp_filtered, 'b');
title('Respiration filtered');
xlabel('Time');
ylabel('z');

subplot(2, 1, 2);
plot(timestamps, heart_filtered, 'r');
title('Heart filtered');
xlabel('Time');
ylabel('X');