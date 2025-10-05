% Ejercicio 1 
clc, clear all, close all 
BR = 32e9;      % Baud Rate
N = 8;          % Oversampling rate
h_taps = 101;   % Pulse shaping taps
fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 consecutive samples at Tx output

% Lista de roll-off a comparar
rolloff_list = [0, 0.1, 0.25, 0.5, 1.0]; 

% FFT params
NFFT = 2048;
f = (-NFFT/2:NFFT/2-1)*fs/NFFT;

colors = lines(length(rolloff_list));

%% Plot en tiempo
figure; hold on; grid on;
for k = 1:length(rolloff_list)
    ro = rolloff_list(k);
    h_rc = raised_cosine(BR/2, fs, ro, h_taps, 0);
    plot(h_rc, 'LineWidth', 1.5, 'Color', colors(k,:));
end
title('Distintos roll-off'); 
xlabel('Samples'); ylabel('Amplitud')
legend(arrayfun(@(r) sprintf('rolloff=%.2f', r), rolloff_list, 'UniformOutput', false))
set(gcf, 'Position', [50 50 600 400],'Color', 'w');

%% Plot en frecuencia
figure; hold on; grid on;
for k = 1:length(rolloff_list)
    ro = rolloff_list(k);
    h_rc = raised_cosine(BR/2, fs, ro, h_taps, 0);
    H_RC = fftshift(abs(fft(h_rc, NFFT)));
    plot(f/1e9, H_RC, 'LineWidth', 1.5, 'Color', colors(k,:));
end
title('Distintos roll-off'); 
xlabel('Freq [GHz]'); ylabel('Amplitud')
legend(arrayfun(@(r) sprintf('rolloff=%.2f', r), rolloff_list, 'UniformOutput', false))
set(gcf, 'Position', [700 50 600 400],'Color', 'w');
